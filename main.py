from astronomy import Equator, Horizon, Time, Body, Observer, Refraction, MoonPhase, SearchMoonQuarter, NextMoonQuarter
import calendar
import datetime
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.font_manager as fm
from backgound import getbackground
from moon_symbol import draw_moon_symbol
from astropy.table import Table
from numpy.lib import recfunctions as rfn
'''
Background part: use sun altitude to color the background
draw a image with shape (365, 24*60), 60min *24h for 365 days
'''
# font
font_path = './fonts/Libertinus-7.040/static/OTF/'
font_files = fm.findSystemFonts(fontpaths=font_path)
for font_file in font_files:
    fm.fontManager.addfont(font_file)
plt.rcParams['font.family'] = 'Libertinus Sans'
plt.rcParams['font.size'] = 13

# observer
north40 = Observer(40, 116, 40)  # Beijing
north24 = Observer(24, 116, 5)  # Amoy
north52 = Observer(52, 4, 20)  # Leiden

obs = north24
time_zone = (obs.longitude) / 15

# year number in Common Era (CE)
year = datetime.date.today().year

# days count in this year
if calendar.isleap(year):
    days_in_year = 366
else:
    days_in_year = 365

# time orgin: noon (12 am) of the first day
time_origin = Time.Make(year, 1, 1, 12 - time_zone, 00, 00.0)

# y scale for all
y_scale_factor = 1 / 14

# font size
fontsize = 8.5

# load lunar mare shape file
lms = []
for f in ['lm1.txt', 'lm2.txt', 'lm3.txt']:
    lm = np.loadtxt('./data/' + f)
    lm[:, 1] = -lm[:, 1]  # counteracting the final y axis flipping
    lms.append(lm)

# read rise and set times, sun alts and EOT
data_path = './data/' + str(year) + '/'
EOT_path = './data/' + str(year) + '/'

rise_table = Table.read(data_path + str(obs.latitude) + '/' + 'rise.fits')
set_table = Table.read(data_path + str(obs.latitude) + '/' + 'set.fits')
transit_table = Table.read(data_path + str(obs.latitude) + '/' +
                           'transit.fits')
twilight_table = Table.read(data_path + str(obs.latitude) + '/' +
                            'twilight.fits')
EOT = np.load(EOT_path + 'EOT.npy')


def draw_moons(ax, moon_rise_or_set):
    moon_rise_or_set_night = moon_rise_or_set[night_mask(moon_rise_or_set)[:,
                                                                           0]]
    for day, hour in moon_rise_or_set_night:
        time_ut = Time.Make(year, 1, 1, 12 - time_zone, 00,
                            00.0).ut + (day - 1) + (hour + 12) / 24
        time = Time(ut=time_ut)
        moon_phase = MoonPhase(time)
        draw_moon_symbol(ax=ax,
                         time=time,
                         moon_phase=moon_phase,
                         lms=lms,
                         day=day,
                         hour=hour,
                         y_zoom_factor=1 / y_scale_factor,
                         zoom_factor=0.075,
                         mares=False)


def draw_patch(ax, vertex, color=[1] * 3, zorder=1):
    codes = np.ones(len(vertex),
                    dtype=mpath.Path.code_type) * mpath.Path.LINETO
    codes[0] = mpath.Path.MOVETO
    path = mpath.Path(vertex, codes)
    patch = mpatches.PathPatch(path,
                               facecolor=color,
                               linewidth=0,
                               zorder=zorder)
    ax.add_patch(patch)


def month_background(sun_rise_or_set_day_hour, left=True):
    if left:
        shift = -1
    else:
        shift = 1
    blue_vertex_1 = np.copy(sun_rise_or_set_day_hour[:, ::-1])
    # day is y coordinate, hour is x coord so here we need to filp it.
    blue_vertex_2 = np.copy(sun_rise_or_set_day_hour[::-1, ::-1])
    blue_vertex_2[:, 0] += shift
    blue_vertex = np.vstack([blue_vertex_1, blue_vertex_2])
    draw_patch(ax,
               blue_vertex,
               color=[226 / 255, 246 / 255, 248 / 255],
               zorder=10)


# get data in astropy table
def rise_set(body: str):
    return [
        rfn.structured_to_unstructured(rise_table['Day', body].as_array()),
        rfn.structured_to_unstructured(set_table['Day', body].as_array())
    ]


def dusk_dawn():
    return [
        rfn.structured_to_unstructured(twilight_table['Day',
                                                      'Dusk'].as_array()),
        rfn.structured_to_unstructured(twilight_table['Day',
                                                      'Dawn'].as_array())
    ]


def transit(body: str):
    return rfn.structured_to_unstructured(transit_table['Day',
                                                        body].as_array())


def night_mask(day_hour, bleed=0.1):

    mask = np.tile(
        np.logical_and(day_hour[:, 1] >= sun_set[:, 1] - bleed,
                       day_hour[:, 1] <= sun_rise[:, 1] + bleed), (2, 1)).T
    return mask


def plot_object_label(ax,
                      body,
                      x,
                      y,
                      verticalalignment='center',
                      horizontalalignment='center',
                      fontsize=fontsize + 2,
                      zorder=100,
                      c='Red'):
    ax.text(x,
            y,
            s=body,
            verticalalignment=verticalalignment,
            horizontalalignment=horizontalalignment,
            fontsize=fontsize,
            zorder=zorder,
            c=c)


def plot_rise_set(ax, bodies: list, zorder, lw, colors):
    for body, c in zip(bodies, colors):
        rise, set = rise_set(body)
        rise = np.ma.array(rise,
                           mask=np.logical_not(night_mask(
                               rise, bleed=0.1))).filled(np.nan)
        set = np.ma.array(set, mask=np.logical_not(night_mask(
            set, bleed=0.1))).filled(np.nan)

        ax.plot(rise[:, 1], rise[:, 0], zorder=zorder, lw=lw, c=c)
        ax.plot(set[:, 1], set[:, 0], zorder=zorder, lw=lw, c=c)


def plot_rise(ax, bodies: list, zorder, lw, colors, object_lable=False):
    for body, c in zip(bodies, colors):
        rise, _ = rise_set(body)
        mask = night_mask(rise, bleed=0.1)
        rise = np.ma.array(rise, mask=np.logical_not(mask)).filled(np.nan)
        ax.plot(rise[:, 1], rise[:, 0], zorder=zorder, lw=lw, c=c)
        if object_lable is True:
            plot_object_label(ax,
                              body,
                              x=rise[mask][1],
                              y=rise[mask][0],
                              verticalalignment='center',
                              horizontalalignment='center',
                              c='Red')
            plot_object_label(ax,
                              body,
                              x=rise[mask][-1],
                              y=rise[mask][-2],
                              verticalalignment='center',
                              horizontalalignment='center',
                              c='Red')


def plot_transit(ax, bodies: list, zorder, lw, colors, object_lable=False):
    for body, c in zip(bodies, colors):
        transit_list = transit(body)
        mask = night_mask(transit_list, bleed=0.1)
        transit_list = np.ma.array(transit_list,
                                   mask=np.logical_not(mask)).filled(np.nan)
        ax.plot(transit_list[:, 1],
                transit_list[:, 0],
                zorder=zorder,
                lw=lw,
                c=c)

        if object_lable is True:
            plot_object_label(ax,
                              body,
                              x=transit_list[mask][1],
                              y=transit_list[mask][0],
                              verticalalignment='center',
                              horizontalalignment='center',
                              c='Red')
            plot_object_label(ax,
                              body,
                              x=transit_list[mask][-1],
                              y=transit_list[mask][-2],
                              verticalalignment='center',
                              horizontalalignment='center',
                              c='Red')


# sun
sun_rise, sun_set = rise_set('Sun')

# for conviniance
sun_rise_days = sun_rise[:, 0]
sun_rise_hours = sun_rise[:, 1]
sun_set_days = sun_set[:, 0]
sun_set_hours = sun_set[:, 1]

# moon
moon_rise, moon_set = rise_set('Moon')

# 绘制
fig, ax = plt.subplots(1, 1, figsize=(15, 15), dpi=300)

# left_white
left_white_vertex = sun_set[:, ::-1]
# day is y coordinate, hour is x coord so here we need to filp it.
left_edge = np.array([[-12, days_in_year], [-12, 1]])
left_white_vertex = np.vstack([left_white_vertex, left_edge])
draw_patch(ax, left_white_vertex, zorder=10)

# right_white
right_white_vertex = sun_rise[:, ::-1]
# day is y coordinate, hour is x coord so here we need to filp it.
right_edge = np.array([[12, days_in_year], [12, 1]])
right_white_vertex = np.vstack([right_white_vertex, right_edge])
draw_patch(ax, right_white_vertex, zorder=10)

# month and date background blue
month_background(sun_set, left=True)
month_background(sun_rise, left=False)

# left and right black lines
black_line_zorder = 20
ax.plot(sun_rise_hours, sun_rise_days, c='k', lw=1, zorder=black_line_zorder)
ax.plot(sun_set_hours, sun_set_days, c='k', lw=1, zorder=black_line_zorder)

ax.plot(sun_rise_hours + 1,
        sun_rise_days,
        c='k',
        lw=0.5,
        zorder=black_line_zorder)
ax.plot(sun_set_hours - 1,
        sun_set_days,
        c='k',
        lw=0.5,
        zorder=black_line_zorder)

# top and butomm line
ax.plot([sun_set_hours[1] - 1, sun_rise_hours[1] + 0.985], [1, 1],
        c='k',
        lw=1,
        zorder=black_line_zorder)
ax.plot([sun_set_hours[-1] - 0.985, sun_rise_hours[-1] + 0.985],
        [days_in_year, days_in_year],
        c='k',
        lw=1,
        zorder=black_line_zorder)

# plot days of date
date_zorder = 30
days_7 = np.arange(1, days_in_year, 7)  # days with step of 7

for day in days_7:
    set_hour = sun_set[day - 1, 1]
    _, _, set_day, _, _, _ = Time.Make(year, 1, day, 12, 00, 00.0).Calendar()
    rise_hour = sun_rise[day - 1, 1]
    _, _, rise_day, _, _, _ = Time.Make(year, 1, day + 1, 12, 00,
                                        00.0).Calendar()

    x_shift = 0.3
    y_shift = 0.5
    ax.text(x=set_hour - x_shift,
            y=day + y_shift,
            s=str(set_day),
            verticalalignment='center',
            horizontalalignment='center',
            fontsize=fontsize,
            zorder=date_zorder)
    ax.text(x=rise_hour + x_shift,
            y=day + y_shift,
            s=str(rise_day),
            verticalalignment='center',
            horizontalalignment='center',
            fontsize=fontsize,
            zorder=date_zorder)

# month line
for month in np.arange(2, 13):
    set_day = int((Time.Make(year, month, 1, 12 - time_zone, 1, 00.0).ut -
                   time_origin.ut) // 1 + 1)
    set_hour = sun_set[set_day - 1, 1]
    ax.plot([set_hour - 1, set_hour], [set_day, set_day],
            lw=0.5,
            c='Black',
            zorder=black_line_zorder)
    rise_day = set_day + 1
    rise_hour = sun_rise[rise_day - 1, 1]
    ax.plot([rise_hour, rise_hour + 1], [rise_day, rise_day],
            lw=0.5,
            c='Black',
            zorder=black_line_zorder)

# sky background
ax.imshow(getbackground(obs, withmoon=False),
          zorder=-1,
          extent=[-12, 12, days_in_year + 0.5, 0 + 0.5])

# 均时差 equation of time
EOT = np.load(data_path + 'EOT.npy')
days_number = np.arange(days_in_year) + 1
ax.plot(EOT, days_number, color='white', lw=0.5, zorder=1)
# EOT axis
ax.plot([0, 0], [1, days_in_year], color='white', lw=0.5, zorder=1)

# grid
grid_c = [0.75] * 3
grid_s = 0.5
min_set_hour = -int(-np.min(sun_set_hours)) - 1
max_set_hour = int(np.max(sun_rise_hours)) + 1

xs = np.arange(min_set_hour, max_set_hour, y_scale_factor)
ys = np.arange(1, days_in_year + 1)

for y in days_7:
    ax.scatter(xs, [y] * xs.shape[0],
               c=[grid_c],
               zorder=0,
               s=grid_s,
               edgecolors='none',
               lw=0)

for x in np.arange(-12, 13, 0.5):
    ax.scatter([x] * ys.shape[0],
               ys,
               c=[grid_c],
               zorder=0,
               s=grid_s,
               edgecolors='none',
               lw=0)

# twilight line
for twilght in dusk_dawn():
    ax.plot(twilght[:, 1],
            twilght[:, 0],
            zorder=1,
            lw=0.5,
            c=[0.75] * 3,
            ls=(0, (8, 5)))

# planet lines
planets_colors_3 = ['firebrick', 'mediumorchid', 'gold']
planets_colors_5 = ['aqua', 'darkorange'] + planets_colors_3
plot_rise_set(ax, ['Mercury', 'Venus', 'Mars', 'Jupiter', 'Saturn'],
              zorder=5,
              lw=0.8,
              colors=planets_colors_5)
plot_transit(ax, ['Mars', 'Jupiter', 'Saturn'],
             zorder=5,
             lw=0.8,
             colors=planets_colors_3)

# stars lines
plot_rise(ax, ['Star1', 'Star2', 'Star3', 'Star4', 'Star5'],
          zorder=5,
          lw=0.5,
          colors=['White'] * 5,
          object_lable=True)
plot_transit(ax, ['Star1', 'Star2', 'Star3', 'Star4', 'Star5'],
             zorder=5,
             lw=0.5,
             colors=['White'] * 5,
             object_lable=True)

# moon phase
draw_moons(ax, moon_rise)
draw_moons(ax, moon_set)

# 朔 望 上下弦
mq = SearchMoonQuarter(time_origin)
quarter_day = (mq.time.ut - time_origin.ut) // 1 + 1.0
quarter_day = int(quarter_day)

moon_rise_masked = moon_rise[np.logical_and(
    moon_rise[:, 1] >= sun_set[:, 1] - 2,
    moon_rise[:, 1] <= sun_rise[:, 1] + 2)]
moon_set_masked = moon_set[np.logical_and(
    moon_set[:, 1] >= sun_set[:, 1] - 2, moon_set[:, 1] <= sun_rise[:, 1] + 2)]

while quarter_day <= days_in_year:
    moon_quarter_rise = moon_rise_masked[moon_rise_masked[:, 0] == quarter_day]
    moon_quarter_set = moon_set_masked[moon_set_masked[:, 0] == quarter_day]
    if moon_quarter_rise.size > 0: # rise 
        draw_moon_symbol(ax=ax,
                        time=mq.time,
                        moon_phase=90 * mq.quarter + 0.0001,
                        lms=lms,
                        day=quarter_day,
                        hour=moon_quarter_rise[0, 1],
                        y_zoom_factor=1 / y_scale_factor,
                        zoom_factor=0.1,
                        mares=True)
    if moon_quarter_set.size > 0: # set
        draw_moon_symbol(ax=ax,
                            time=mq.time,
                            moon_phase=90 * mq.quarter + 0.0001,
                            lms=lms,
                            day=quarter_day,
                            hour=moon_quarter_set[0, 1],# set
                            y_zoom_factor=1 / y_scale_factor,
                            zoom_factor=0.1,
                            mares=True)

    mq = NextMoonQuarter(mq)
    quarter_day = (mq.time.ut - time_origin.ut) // 1 + 1

# fig basic adjustment
ax.set_xlim(-12, 12)
ax.set_ylim(days_in_year, 1)

ax.axis('off')

# output scale
ax.set_aspect(y_scale_factor)

# plt.show()
fig.savefig('test_' + str(obs.latitude) + '.pdf')
