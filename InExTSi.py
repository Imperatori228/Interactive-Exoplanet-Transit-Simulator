import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import batman

# ვარსკვლავის სისტემის განმსაზღვრელი პარამეტრები და ფიზიკური კონსტანტები - Constants and Star System Definitions
R_SUN = 6.957e8; M_SUN = 1.989e30; R_EARTH = 6.371e6; AU = 1.496e11
STAR_SYSTEMS = {
    "Sun-like (G-type)": {
        "star_radius_m": R_SUN, "star_mass_kg": M_SUN, "limb_dark_u": [0.4, 0.2], "color": "gold",
        "hz_inner_au": 0.95, "hz_outer_au": 1.67
    },
    "M-dwarf": {
        "star_radius_m": 0.12 * R_SUN, "star_mass_kg": 0.08 * M_SUN, "limb_dark_u": [0.6, 0.3], "color": "red",
        "hz_inner_au": 0.03, "hz_outer_au": 0.06
    }
}
# გამოთვლების აჩქარებისთვის ცვლილების წერტილების რაოდენობის შემცირება - Reduce the number of points for faster calculation
NUM_TIME_POINTS = 300
T_HOURS = np.linspace(-5, 5, NUM_TIME_POINTS)

def calculate_transit_signal(star_params, orbital_distance_m, planet_radius_m, atmosphere_height_m):
    
    star_radius_m = star_params["star_radius_m"]
    star_mass_kg = star_params["star_mass_kg"]
    G = 6.67430e-11
    P_seconds = np.sqrt(4 * np.pi**2 * orbital_distance_m**3 / (G * star_mass_kg))
    P_days = P_seconds / (60 * 60 * 24)
    params = batman.TransitParams()
    params.t0 = 0.; params.per = P_days; params.a = orbital_distance_m / star_radius_m
    params.inc = 90.; params.ecc = 0.; params.w = 90.
    params.limb_dark = "quadratic"; params.u = star_params["limb_dark_u"]
    t_days = T_HOURS / 24.0
    params.rp = planet_radius_m / star_radius_m
    model_solid = batman.TransitModel(params, t_days); flux_solid = model_solid.light_curve(params)
    depth_solid_ppm = (1 - min(flux_solid)) * 1e6
    planet_radius_with_atmos = planet_radius_m + atmosphere_height_m
    params.rp = planet_radius_with_atmos / star_radius_m
    model_atmos = batman.TransitModel(params, t_days); flux_atmos = model_atmos.light_curve(params)
    depth_atmos_ppm = (1 - min(flux_atmos)) * 1e6
    atmospheric_signal_ppm = depth_atmos_ppm - depth_solid_ppm
    return flux_atmos, atmospheric_signal_ppm, orbital_distance_m / AU

# პლოტის შექმნა - Setup the plot
fig, ax = plt.subplots(figsize=(10, 8))
plt.subplots_adjust(bottom=0.4)
ax.set_title('სიმულატორი (Simulator): M-კლასის ჯუჯა ვარსკვლავის უპირატესობა (The M-Dwarf Advantage)', fontsize=16)
ax.set_xlabel('დრო ტრანზიტის ცენტრიდან (საათებში) - Time from Transit Center (Hours)'); ax.set_ylabel('Relative Stellar Flux')
ax.grid(True, linestyle='--', alpha=0.6)

# პლოტის ხაზების შექმნა, შემდეგ მათი განახლბა უფრო კარგი მუშაობისთვის - Create plot lines
line_sun, = ax.plot(T_HOURS, np.ones(NUM_TIME_POINTS), color=STAR_SYSTEMS["Sun-like (G-type)"]['color'], lw=2, label='მზისნაირი სისტემა - Sun-like System')
line_mdwarf, = ax.plot(T_HOURS, np.ones(NUM_TIME_POINTS), color=STAR_SYSTEMS["M-dwarf"]['color'], lw=2, label='M-ჯუჯა სისტემა - M-dwarf System')
ax.legend(loc='lower center')
results_text_obj = ax.text(0.02, 0.5, "", transform=ax.transAxes, fontsize=10,
                           verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

# სლაიდერების შექმნა განსაზღვრული ნაბიჯებით რათა უფრო კარგად იმუშაოს - Create sliders with discrete steps for faster updates
ax_planet_radius = plt.axes([0.25, 0.22, 0.65, 0.03])
slider_planet_radius = Slider(ax=ax_planet_radius, label='პლანეტის რადიუსი (Planet Radius) [R⊕]', valmin=0.5, valmax=2.5, valinit=1.0, color='cyan', valstep=0.05)

ax_atmos_height = plt.axes([0.25, 0.15, 0.65, 0.03])
slider_atmos_height = Slider(ax=ax_atmos_height, label='ატმოს. სიმაღლე (Atmos. Height) [km]', valmin=10, valmax=200, valinit=50, color='lightblue', valstep=5)

ax_hz_pos = plt.axes([0.25, 0.08, 0.65, 0.03])
slider_hz_pos = Slider(ax=ax_hz_pos, label='ცხოვრებისეული ზონის პოზიცია (HZ Position)', valmin=0, valmax=1, valinit=0.5, color='lightgreen', valstep=0.05)

def update_plot(val):
    # სლაიდერებიდან მნიშვნელობების მიღება - Get values from sliders
    planet_radius_re = slider_planet_radius.val
    atmosphere_height_km = slider_atmos_height.val
    hz_pos = slider_hz_pos.val

    # სი სისტემის ერთეულებში გადაყვანა - Convert to SI units
    planet_radius_m = planet_radius_re * R_EARTH
    atmosphere_height_m = atmosphere_height_km * 1000

    # ორბიტალური დისტანციის გამოთვლა - Calculate orbital distances
    sun_params = STAR_SYSTEMS["Sun-like (G-type)"]
    sun_dist_au = sun_params["hz_inner_au"] + hz_pos * (sun_params["hz_outer_au"] - sun_params["hz_inner_au"])
    sun_dist_m = sun_dist_au * AU

    mdwarf_params = STAR_SYSTEMS["M-dwarf"]
    mdwarf_dist_au = mdwarf_params["hz_inner_au"] + hz_pos * (mdwarf_params["hz_outer_au"] - mdwarf_params["hz_inner_au"])
    mdwarf_dist_m = mdwarf_dist_au * AU

    # სიმულაციას რთავს - Run simulations
    sun_flux, sun_signal, _ = calculate_transit_signal(sun_params, sun_dist_m, planet_radius_m, atmosphere_height_m)
    mdwarf_flux, mdwarf_signal, _ = calculate_transit_signal(mdwarf_params, mdwarf_dist_m, planet_radius_m, atmosphere_height_m)

    # პლოტის მონაცემების განახლება - Update the plot data
    line_sun.set_ydata(sun_flux)
    line_mdwarf.set_ydata(mdwarf_flux)
    ax.relim() # ღერძის ზღვარის გადათვლა - Recalculate axis limits
    ax.autoscale_view(scalex=False, scaley=True) # მხოლოდ Y ღერძის ზომის შეცვლა - Rescale Y axis only

    # რეზულტატების ტექსტის განახლება - Update results text
    advantage_factor = mdwarf_signal / sun_signal if sun_signal > 0 else float('inf')
    results_text = (
        f"--- სიმულაციის რეზულტატები (Simulation Results) ---\n"
        f"პლანეტის რადიუსი (Planet Radius): {planet_radius_re:.2f} R⊕\n"
        f"ატმოს. სიმაღლე (Atmos. Height): {atmosphere_height_km:.0f} km\n"
        f"ცხოვრებისეული ზონის პოზიცია (HZ Position): {'Inner' if hz_pos < 0.33 else 'Mid' if hz_pos < 0.66 else 'Outer'}\n\n"
        f"ორბიტალური დისტანცია (Orbital Dist) (მზე - Sun): {sun_dist_au:.2f} AU\n"
        f"ორბიტალური დისტანცია (Orbital Dist) (M-ჯუჯა - M-dwarf): {mdwarf_dist_au:.3f} AU\n\n"
        f"ატმოს. სიგნალი (Atmos. Signal) (Sun-like): {sun_signal:.2f} ppm\n"
        f"ატმოს. სიგნალი (Atmos. Signal) (M-dwarf): {mdwarf_signal:.2f} ppm\n\n"
        f"უპირატესობის ფაქტორი (Advantage Factor): {advantage_factor:.1f}x"
    )
    results_text_obj.set_text(results_text)
    
    fig.canvas.draw_idle()

# სლაიდერების დაკავშირება და ჩართვა - Link sliders and run
slider_planet_radius.on_changed(update_plot)
slider_atmos_height.on_changed(update_plot)
slider_hz_pos.on_changed(update_plot)
update_plot(None)

# პროგრამის მთელს ეკრანზე გადიდებისთვის - To maximize the window
try:
    manager = plt.get_current_fig_manager()
    if hasattr(manager, 'window'):
        # For Qt backend
        if hasattr(manager.window, 'showMaximized'):
            manager.window.showMaximized()
        # For Tk backend
        elif hasattr(manager.window, 'state'):
            manager.window.state('zoomed')
except (AttributeError, ImportError):
    print("Could not automatically maximize the plot window. You may need to do it manually.")
plt.show()
