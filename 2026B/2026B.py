# # Computational Mathematics
# ## Project 2026B
# ## Candidate Number: 1105769

# General references and documentation used throughout:
# - https://docs.python.org/3.12/library/stdtypes.html
# - https://numpy.org/doc/stable/reference/index.html
# - https://matplotlib.org/stable/api/index.html

# Imports for the script

from astropy.time import Time, TimeDelta

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt

import numpy as np
import sympy as sp

# To render sympy equations
from publish import render

# Just to suppress the warnings so the report looks cleaner
# https://docs.python.org/3/library/warnings.html#warnings.filterwarnings
from cartopy.io import DownloadWarning
import warnings
warnings.filterwarnings("ignore", category=DownloadWarning)

# ***
# ### Question B.1.
#
# > Use the `numpy.polynomial.Polynomial` class to represent the
# > different Besselian elements. Plot the coordinates
# > $(x(t), y(t))$ on the Besselian plane for the elements given in Table B.1,
# > for $(t - t_0) \in [-2.5, 2.5]$ hours.
# >
# > [1 mark]

t0 = Time("1715-05-03 10:00:00", scale="tt")
delta_T = 9.6 # seconds

# $(x, y)$ is the centre of the shadow on the fundamental plane
x = np.polynomial.Polynomial([0.0718130, 0.5682352, 0.0000137, -0.0000096])
y = np.polynomial.Polynomial([0.7433290, 0.1231166, -0.0001459, -0.0000020])

# $l_1$ is the radius of the penumbral shadow
l1 = np.polynomial.Polynomial([0.5330720, 0.0000319, -0.0000128])

# $l_2$ is the radius of the umbral shadow
l2 = np.polynomial.Polynomial([-0.0130010, 0.0000318, -0.0000127])

# $\d$ is the declination of the shadow axis
d = np.pi / 180 * np.polynomial.Polynomial([15.5370197, 0.0120410, -0.0000030]) # (radians)

# $\mu$ is the hour angle
mu = np.pi / 180 * np.polynomial.Polynomial([330.845947, 15.002640]) # (radians)

# $f_1$ is the semi-vertex angle of the penumbral cone
# $f_2$ is the semi-vertex angle of the umbral cone
# They change so slowly that they are treated as a constant.
tan_f1 = 0.0046323
tan_f2 = 0.0046092

def question_b1():
    fig, ax = plt.subplots(layout="constrained")

    # Generate a sample of (t - t_0) to evaluate x and y at
    ts = np.linspace(-2.5, 2.5, 100_000)

    # Calculate x and y at those times
    xs = x(ts)
    ys = y(ts)

    # Plot the coordinates
    ax.plot(xs, ys, label="$(x(t - t_0), y(t - t_0))$", color="C0")

    # Mark where was t - t_0 = -2.5 and t - t_0 = 2.5 on the graph to give
    # a sense of direction.
    # Use zorder to make the dots appear above the lines
    ax.scatter([xs[0], xs[-1]], [ys[0], ys[-1]], color="C1", zorder=2)
    ax.annotate("$t - t_0 = -2.5$", (xs[0], ys[0]), ha="right", textcoords="offset points", xytext=(-5, 0))
    ax.annotate("$t - t_0 = 2.5$", (xs[-1], ys[-1]), textcoords="offset points", xytext=(5, 0))

    # Just adding chart elements for the graph
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(r"Centre of the shadow on the fundamental plane for $(t - t_0) \in [-2.5, 2.5]$")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal") # Set equal aspect ratio for x and y axes

    # Extend axis limits to allow annotations to fit better
    # https://stackoverflow.com/a/32693077
    ax.margins(x=0.3, y=0.25)

    fig.show()

question_b1()

# The graph looks very linear, which absolutely makes sense as the $t^2$
# coefficients are very small compared to the $t$ coefficients
# ($0.5682352 \gg 0.0000137$ and $0.1231166 \gg -0.0001459$).


# ***
# ### Question B.2.
#
# > Plot the difference $\theta - \phi$ as a function of geocentric
# > latitude $\theta$. On the Earth, what is the maximal value of $|\theta - \phi|$, and
# > where does it occur? Report your results in both radians and degrees.
# > What is the maximal value of the difference on Jupiter, where $a \approx 71492 km$
# > and $b \approx 66854 km$? Comment on your results.
# >
# > [2 marks]

a_earth = 6378137.0 # (m)
b_earth = 6356752.314245 # (m)
a_jupiter = 71492e3 # (m)
b_jupiter = 66854e3 # (m)

def question_b2():
    # These are the values of theta (geocentric latitude) which we will use to
    # calculate the difference between theta and phi (geographic latitude).
    # I decided on 100_000 points which should be more than sufficient to
    # make the graph ``look'' continuous.
    # Note I kept thetas in radians as the formula for computing phi using theta
    # is easier to compute from radians
    thetas = np.linspace(-np.pi/2, np.pi/2, 100_000) # (radians)

    # Setup graph area to allow two plots side by side
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")


    # Using a loop as we are plotting essentially the same features for both
    # the Earth and Jupiter
    planet_radii = {
        "Earth": (a_earth, b_earth),
        "Jupiter": (a_jupiter, b_jupiter)
    }
    # We also need the index to know which plot to be plotting in,
    # so we use enumerate
    for idx, (planet, (a, b)) in enumerate(planet_radii.items()):
        # Compute the geographic latitude for each value of geocentric latitude
        phi = np.arctan((a / b)**2 * np.tan(thetas)) # (radians) given in (B.2.3b)

        # Using np.degrees to convert radians to degrees so that the code
        # is slightly more readable.
        # I chose to plot in degrees as that seems to be the usual unit used
        # for latitude and longitude diagrams.
        axs[idx].plot(np.degrees(thetas), np.degrees(thetas - phi))

        # Add chart elements
        axs[idx].set_xlabel(r"$\theta ({}^\circ)$")
        axs[idx].set_ylabel(r"$\theta − \phi ({}^\circ)$")
        axs[idx].set_title(f"Difference between Geocentric ($\\theta$) and \n Geographic ($\\phi$) latitude against $\\theta$ for {planet}")
        axs[idx].grid(True, alpha=0.3)
        # Setting xticks at multiples of 22.5 degrees rather than the default 25
        # so that it better shows how close the maximum occurs around 45 degrees
        axs[idx].set_xticks(np.linspace(-90, 90, 9))

        # Find the index with the maximal difference
        max_diff_idx = np.argmax(np.abs(thetas - phi))
        # Get the corresponding difference at that index
        max_diff_deg = np.degrees(thetas[max_diff_idx] - phi[max_diff_idx]) # (degrees)
        # The value of theta that results in the maximal difference
        theta_max_deg = np.degrees(thetas[max_diff_idx]) # (degrees)

        # Print the maximal difference and when it occurs
        print(f"On {planet}, max(|theta - phi|) = {np.radians(max_diff_deg): .5f} rad = {max_diff_deg: .3f} deg")
        print(f"at theta = {thetas[max_diff_idx]: .5f} rad = {theta_max_deg: .3f} deg")

        # Mark the point of maximum absolute difference on the graph.
        # Setting zorder=2 so that this point appears "in front" of
        # the line plot rather than "behind"
        axs[idx].scatter(theta_max_deg, max_diff_deg, color="C1", zorder=2)
        # Label that point on the graph
        axs[idx].annotate(fr"$\max|\theta - \phi| = {max_diff_deg: .3f}^\circ$ at $\theta = {theta_max_deg: .2f}^{{\circ}}$",
                          xy=(theta_max_deg, max_diff_deg),
                          xytext=(0, 7), # shift the label slightly up so that it is less cramped
                          textcoords="offset points")

        # I want to draw a dotted line from this point of maximal difference to
        # the x-axis, but I want to keep the lower y limit the same, without the
        # plot area extending further down to accomodate the extra line I am
        # drawing for visual purposes.
        # I also want to increase the upper y limit so that the label doesn't
        # "look" cramped in the plot
        ymin_plot, ymax_plot = axs[idx].get_ylim()
        axs[idx].set_ylim(ymin_plot, ymax_plot * 1.05)
        axs[idx].vlines(theta_max_deg, ymin=ymin_plot, ymax=max_diff_deg, colors="C1", linestyles="dotted", alpha=0.7)

    plt.show()

question_b2()

# We can note that $\theta - \phi$ is zero at the equator and the poles.
# This maximum occurs near $\frac{\pi}{4} \text{rad} = 45^\circ$, but slightly closer to the equator than the pole.
#
# This agrees with some sources on the internet, e.g.
# - https://www.johndcook.com/blog/2011/09/17/latitude-doesnt-exactly-mean-what-i-thought/
# > The maximum occurs at 44.9 degrees and equals 0.1917.
# 
# Even though we found a negative $\theta$, but our function is odd wrt $\theta$, hence a maximal difference at $\theta$ would also imply that it is maximal at $-\theta$.
# It just so happened, due to numerical precision, that argmax found the negative value to have a larger difference.
#
# We also found that the maximum difference is much greater for Jupiter ($0.06$ rad) compared to the Earth ($0.003$ rad).
# - https://www.universetoday.com/articles/jupiter-compared-to-earth
# > Also like Earth, Jupiter's shape is that of an oblate spheroid. In fact, Jupiter's polar flattening is greater than that of Earth's - 0.06487 ± 0.00015 compared to 0.00335.
#
# These numbers do agree with the ones we computed here, and the reasoning for this phenomenon is explained by:
# - https://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node52.html
# > This degree of flattening is much larger than that of the Earth, owing to Jupiter's relatively large radius (about ten times that of Earth), combined with its relatively short rotation period (about 0.4 days).


# ***
# ### Question B.3.
#
# > Write a function `compute_latitude_longitude` that implements the projection
# > of points $(\sigma, \eta)$ on the fundamental plane along the shadow axis
# > given by (B.3.2)–(B.3.7). For later use, the function should take $\sigma$
# > and $\eta$ as matrices of shape $n_{pts} \times n_{times}$, where $n_{pts}$
# > is the number of points at a given time $t - t_0$ and $n_{times}$ is the
# > number of different times to consider.
# >
# > What is the geographical latitude and longitude of the centre line at
# > the time of maximal eclipse, $t_{max} =$ 1715-05-03 09:39:30 TDT?
# >
# > [4 marks]
#
# > $$\tilde\rho = \sqrt{1 - e^2 \cos^2 d}, \sin\tilde{d} = \frac{\sin d}{\tilde\rho}, \cos\tilde{d} = \frac{(1 - f)\cos d}{\tilde\rho}, \tag{B.3.1}$$
# > $$\tilde\eta = \frac{\eta}{\tilde\rho}. \tag{B.3.2}$$
# > $$\tilde\zeta = \sqrt{1 - \sigma^2 - \tilde\eta^2}, \text{for } 1 - \sigma^2 - \tilde{\eta}^2 \ge 0. \tag{B.3.3}$$
# > $$\sigma(x)^2 - \tilde{\eta}(y)^2 = 1. \tag{B.3.4}$$
# > $$\begin{bmatrix}\cos\tilde\phi \sin\theta \\ \sin\tilde\phi \\ \cos\tilde\phi \cos\theta \end{bmatrix} = \begin{bmatrix} 1 & 0 & 0 \\ 0 & \cos\tilde{d} & \sin\tilde{d} \\ 0 & -\sin\tilde{d} & \cos\tilde{d} \end{bmatrix}\begin{bmatrix}\sigma \\ \tilde\eta \\ \tilde\zeta \end{bmatrix}. \tag{B.3.5}$$
# > $$\lambda = \theta - \mu - \Omega\Delta T. \tag{B.3.6}$$
# > $$\tan\phi = \frac{\tan\tilde\phi}{1 - f}. \tag{B.3.7}$$

# some constants that appear in the calculatione below:
f_earth = 1 - b_earth / a_earth # (B.2.1) flattening
e2_earth = 2 * f_earth - f_earth**2 # (B.2.2) eccentricity

# Constant encoding the difference between TDT and UT
Omega = 7.29211585e-5  # rad/s

def compute_latitude_longitude(sigma, eta, ts):
    # Using assert statements to ensure function input is in right format
    assert sigma.shape == eta.shape, "sigma and eta must have the same shape (npts, ntimes)!"
    assert ts.shape[0] == sigma.shape[1], "sigma must have shape (npts, ntimes), and ts must have shape (ntimes,)!"

    # Calculating auxiliary values by (B.3.1)
    ds = d(ts)
    rho_tilde = np.sqrt(1 - e2_earth * np.cos(ds)**2)
    sin_d_tilde = np.sin(ds) / rho_tilde
    cos_d_tilde = (1 - f_earth) * np.cos(ds) / rho_tilde

    # Normalising y-coordinate by (B.3.2)
    eta_tilde = eta / rho_tilde

    # Calculate auxiliary depth by (B.3.3)
    # If $1 - \tilde\sigma^2 - \tilde\eta^2 < 0$, `np.sqrt` will return `np.nan`
    zeta_tilde = np.sqrt(1 - sigma**2 - eta_tilde**2)

    # Transform cartesian coordinates into observer's normalised latitude and local hour angle by (B.3.5)
    # Implementing this as a matrix multiplication is more troublesome,
    # but we can simply expanded out
    cos_phi_tilde_sin_theta = sigma
    sin_phi_tilde = cos_d_tilde * eta_tilde + sin_d_tilde * zeta_tilde
    cos_phi_tilde_cos_theta = -sin_d_tilde * eta_tilde + cos_d_tilde * zeta_tilde

    # Solving for theta (geocentric latitude)
    theta = np.arctan2(cos_phi_tilde_sin_theta, cos_phi_tilde_cos_theta)

    # Solving for lambda (longitude)
    lambda_ = theta - mu(ts) - Omega * delta_T # (B.3.6) longitude
    lambda_ = (lambda_ + np.pi) % (2.0 * np.pi) - np.pi # to ensure lambda lies between [-pi, pi]

    # Solving for phi (geographic lattitide)
    cos_phi_tilde = np.sqrt(cos_phi_tilde_sin_theta**2 + cos_phi_tilde_cos_theta**2)
    phi_tilde = np.arctan2(sin_phi_tilde, cos_phi_tilde) # parametric latitude
    phi = np.arctan(np.tan(phi_tilde) / (1 - f_earth)) # geographic latitude

    # returns latitude and longitude in degrees
    return np.degrees(phi), np.degrees(lambda_)


def question_b3():
    # Using t_max given by the question:
    t_max = Time("1715-05-03 09:39:30", scale="tt")
    # Calculate time difference to t_0
    t_diff = (t_max - t0).to_value("hour")

    # Centre line: sigma = x(t), eta = y(t)
    # We are only interested in the latitude and longitude of the single point (centre line) (so npts = 1)
    # at this specific time t_max (so ntimes = 1).
    sigma_max = np.array([[x(t_diff)]])
    eta_max = np.array([[y(t_diff)]])
    ts = np.array([t_diff])
    latitudes, longitudes = compute_latitude_longitude(sigma_max, eta_max, ts)

    print(f"At time {t_max} TDT (maximal eclipse time given by the question):")
    print(f"latitude = {latitudes[0, 0]:.3f} deg, longitude = {longitudes[0, 0]:.3f} deg")

question_b3()

# From our above calculations, the coordinates of the centre line at time 
# of maximal eclipse would be 60.0416°N, 20.0961°E.
# However, this seems to differ from the coordinates given by NASA at
# https://eclipse.gsfc.nasa.gov/SEsearch/SEdata.php?Ecl=+17150503
# > Latitude: 59.4° N, Longitude: 17.9° E
#
# On closer inspection, NASA claims that the maximal eclipse occurs at 09:3**6**:30 TDT,
# as opposed to 09:3**9**:30 TDT as what our question states.
#
# We can perform some analysis to see which might be more accurate.
# According to https://eclipse.gsfc.nasa.gov/SEmono/reference/map.html:
# > Greatest eclipse is defined as the instant when the axis of the Moon's shadow passes closest to Earth's center. 
#
# Hence we would want to find the time $t$ that minimises $x(t)^2 + y(t)^2$.

def find_maximal_eclipse_time():
    # Find t_max by minimising x(t)^2 + y(t)^2 analytically
    # Since x and y are just polynomials, we can find the extrema by finding
    # the roots of the first derivative    
    deriv = (x * x + y * y).deriv()
    roots = deriv.roots()
    print(f"The local extrema of x(t)^2 + y(t)^2 are when t in {roots}")
    print()
    print(f"We can note that out of these roots, the only time for maximal eclipse that is reasonable is {roots[2]:.3f} hours from t_0")
    print()
    t_max_new = roots[2]
    second_deriv = deriv.deriv()
    print(f"d^2/dt^2 (x(t)^2 + y(t)^2) at t = {t_max_new:.3f} is {second_deriv(t_max_new):.3f} > 0")
    print(f"Therefore {t_max_new:.3f} is a minimum indeed by second derivative test")
    print()
    print(f"Maximal eclipse should occur {t_max_new:.3f} hours from t_0 = {t0 + TimeDelta(t_max_new * 3600, format="sec")} TDT")

    sigma_max = np.array([[x(t_max_new)]])
    eta_max = np.array([[y(t_max_new)]])
    ts = np.array([t_max_new])
    latitudes, longitudes = compute_latitude_longitude(sigma_max, eta_max, ts)

    print(f"At time {t0 + TimeDelta(t_max_new * 3600, format="sec")} TDT (maximal eclipse time we just found):")
    print(f"latitude = {latitudes[0, 0]:.3f} deg, longitude = {longitudes[0, 0]:.3f} deg")

find_maximal_eclipse_time()

# Hence the maximal eclipse time and coordinates provided by NASA are likely more accurate.
# Perhaps this was a simple typo in the question.


# ***
# ### Question B.4.
#
# > Using numerical rootfinding, determine the start and
# > end times for the 1715 eclipse. Express your answer in hours from $t_0$.
# >
# > [3 marks]
# > 
# > [Hint: an eclipse on Earth cannot last more than 8 hours. Use this to
# > generate good initial guesses, then apply a standard rootfinding algorithm.]

def shadow_hit_earth(t):
    # 1 - sigma^2 - eta_tilde^2 is positive when the shadow hits the earth,
    # negative if it does not.
    # We use this as an indicator function for the bisection search
    # which should find the time t that the shadow just starts to or stops
    # hitting the earth.

    # We can check that this function is indeed continuous, hence bisection search
    # will work due to Intermediate Value Theorem, given start and end guesses that
    # have opposite signs.
    rho_tilde = np.sqrt(1 - e2_earth * np.cos(d(t))**2)
    eta_tilde = y(t) / rho_tilde
    return 1 - x(t)**2 - eta_tilde**2

def bisection_search(f, low, high, eps=1e-9):
    # Given a continuous function `f`, and two initial guesses `low` and `high` 
    # with opposite signs, find a root using bisection, accurate to uncertainty
    # `eps` from the true root.
    # This is a method with linear convergence rate.

    low_val = f(low)
    high_val = f(high)
    assert low_val * high_val < 0, "f(low) and f(high) must have opposite signs!"

    if low_val > high_val:
        # Swap the lower and upper search bounds so that f(low) \le f(high)
        # Just makes the logic slightly cleaner later
        return bisection_search(f, high, low, eps)

    # From here on WLOG assume that f(high_val) > f(low_val)

    # After $n$ iterations, the length of the interval is $|high - low| / 2^n$.
    # Hence the loop will terminate in at most $n = ceil(log_2(|high - low| / eps))$
    while np.abs(high - low) > eps:
        mid = (low + high) / 2
        mid_val = f(mid)

        if mid_val >= 0:
            high = mid
            high_val = mid_val
        else:
            low = mid
            low_val = mid_val

    return (low + high) / 2

def question_b4():
    # Using the hint that an Eclipse on Earth cannot last more than 8 hours,
    # we are going to guess that at t = -8 the eclipse has not started,
    # at t = 0, the eclipse is ongoing,
    # and at t = 8 the eclipse has ended.

    start_low = -8
    start_high = 0
    end_low = 0
    end_high = 8

    # Times found by bisection method
    start_time = bisection_search(shadow_hit_earth, start_low, start_high)
    end_time = bisection_search(shadow_hit_earth, end_low, end_high)

    # Print the results
    print(f"Start time of eclipse is {start_time :.3f}, end time is {end_time :.3f} (hours from t_0)")
    return start_time, end_time

# Store the start and end times of the eclipse in a variable
# since I need to use it for Question B5.
start_time, end_time = question_b4()


# ***
# ### Question B.5.
#
# > Plot the centre line of the 1715 eclipse on a map of the
# > Earth. The map should be drawn with the `cartopy` package and
# > should render coastline boundaries. Consult the `cartopy` documentation
# > and examples online to learn how to use it. Choose a suitable
# > (i) map projection (ii) colour scheme for land/ocean (iii) resolution
# > the coastline (iv) map region so that the entire centre line is visible.
# >
# > [5 marks]

# References used
# - https://cartopy.readthedocs.io/stable/reference/index.html
# - https://cartopy.readthedocs.io/stable/matplotlib/intro.html
# - https://cartopy.readthedocs.io/stable/gallery/lines_and_polygons/features.html
# - https://cartopy.readthedocs.io/stable/gallery/lines_and_polygons/hurricane_katrina.html
# - https://stackoverflow.com/a/51134324

def get_centre_line_coords(t_start, t_end, n=10_000):
    # Compute centre line latitude and longitudes, and convert into degrees.
    ts = np.linspace(t_start, t_end, n)
    xs = x(ts).reshape(1, -1)
    ys = y(ts).reshape(1, -1)

    latitudes, longitudes = compute_latitude_longitude(xs, ys, ts)
    return latitudes.squeeze(), longitudes.squeeze()

def question_b5():
    # Compute the centre line
    latitudes, longitudes = get_centre_line_coords(start_time, end_time)

    # Drawing the map
    fig = plt.figure(figsize=(12, 7), layout="constrained")
    # I chose Mercator projection, which is quite standard across many uses,
    # including Google Maps.
    # It's also a cylindrical projection, so it's just easier to look at
    # compared to a conic projection (LambertConformal or Albers), and can 
    # show a larger area than a Azimuthal projection (like Orthographic).
    # Furthermore, Edmond Halley draw his map with a cylindrical projection.
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    
    # Show Land, Ocean and Coastline borders with 50m resolution, which seems
    # to be sufficient at this map scale.
    # Also simply used the default color scheme for Land and Oceans
    ax.add_feature(cfeature.LAND.with_scale("50m"))
    ax.add_feature(cfeature.OCEAN.with_scale("50m"))
    ax.add_feature(cfeature.COASTLINE.with_scale("50m"))
    # Show most of the world map (excluding some extreme lattitudes as Mercator
    # tends to inflate them too much, and is also not very important to our map
    # as the centre line doesn't extend that far.)
    ax.set_extent([0, 360, -60, 75], crs=ccrs.PlateCarree())

    # Draw the centre line on the map
    ax.plot(longitudes, latitudes, transform=ccrs.Geodetic(), color="red", label="Centre line")

    # Add chart elements (gridlines, legend, title)
    ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    plt.legend()
    plt.title("1715 Solar Eclipse Centre Line")

question_b5()

# The map does seem about right, comparing with https://photoephemeris.com/en/eclipses/solar/TSE1715


# ***
# ### Question B.6.
#
# > Write a function that compute $N$ points uniformly distributed in angle
# > $\alpha$ on the boundary of the umbral cone in Besselian
# > coordinates at a given time $t - t_0$ using the fixed-point iteration
# > (B.4.1) - (B.4.5). Terminate each fixed point iteration when the update
# > to $\zeta$ is sufficiently small. What is the Besselian coordinates of the
# > point on the umbral boundary with $\alpha = \pi/2$ at $t - t_0 = -0.87$?
# >
# > [2 marks]
#
# > $$L_2 = l_2 - \zeta \tan f_2 \tag{B.4.1}$$
# > $$(\sigma - x(t))^2 + (\eta - y(t))^2 < L_2(t, \zeta)^2 \tag{B.4.2}$$
# > $$\sigma(\alpha, t) = x(t) + L_2(t) \cos\alpha \tag{B.4.3}$$
# > $$\eta(\alpha, t) = y(t) + L_2(t) \sin\alpha \tag{B.4.4}$$
# > $$\sigma^2 + \eta^2 + \frac{\zeta^2}{1 - e^2} = 1 \implies \zeta = \sqrt{(1 - e^2)(1 - \sigma^2 - \eta^2)} \tag{B.4.5}$$

def compute_umbral_boundary(N, t, eps=1e-9):
    # computes N equally angle distant points and finds zeta for each of them
    # using fixed-point iteration.

    # Generate uniformly distributed angles
    alphas = np.linspace(0, 2 * np.pi, N, endpoint=False)

    # Our initial guess is where $(\sigma, \eta) = (x, y)$, which is the position of the centre of the umbral cone.
    # This initial guess is always valid as $x^2 + y^2 < 1$ whenever the shadow axis passes through Earth, i.e. during the duration of the eclipse.
    zetas = np.full_like(alphas, np.sqrt((1 - e2_earth) * (1 - x(t)**2 - y(t)**2)))

    # Performs fixed point iteration for each angle alpha:
    iters = 0
    while True:
        iters += 1
        # Calculate next guess for `zetas` as `zetas_`
        L2 = l2(t) - zetas * tan_f2 # (B.4.1)
        sigma = x(t) + L2 * np.cos(alphas) # (B.4.3)
        eta = y(t) + L2 * np.sin(alphas) # (B.4.4)
        zetas_ = np.sqrt((1 - e2_earth) * (1 - sigma**2 - eta**2)) # (B.4.5)

        # Checks if the update to zeta is less than eps for all alpha
        if np.abs(zetas_ - zetas).max() < eps:
            # If update is small return the next guess as the final zeta
            return alphas, zetas_, iters

        # Else update zeta and repeat fixed-point iteration
        zetas = zetas_

# To see why fixed-point iteration works:
# Note that the function we are looking for the fixed point on is
# $$g(\zeta) = \sqrt{(1 - e^2)(1 - \sigma^2 - \eta^2)}.$$
# We then need to express it completely in terms of $\zeta$ and check its derivative.

def fixed_point_analysis():
    # Define the symbols
    x, y, l2, alpha = sp.symbols(r"x y l_2 \alpha", real=True)
    zeta, e2, tan_f2_ = sp.symbols(r"\zeta e^2 \tan(f_{2})", nonnegative=True)
    # Here we define L2, sigma and eta as expressions of zeta, e^2, x, y, l2, tan f_2, alpha
    L2 = l2 - zeta * tan_f2_ # (B.4.1)
    sigma = x + L2 * sp.cos(alpha) # (B.4.3)
    eta = y + L2 * sp.sin(alpha) # (B.4.4)

    # Our fixed point function
    print("The fixed-point function is:")
    g = sp.sqrt((1 - e2) * (1 - sigma**2 - eta**2)) # (B.4.5)
    render(g, name=r"g(\zeta)")

    # Compute the derivative g'(zeta)
    print("The derivative of g is:")
    g_diff = sp.diff(g, zeta).simplify()
    render(g_diff, name=r"g'(\zeta)")
    # It looks really ugly right now, but we can simplify it through some substitutions.
    sigma_sym, eta_sym, L2_sym = sp.symbols(r"\sigma \eta L_2", real=True)
    # Here we do the reverse substitution of the expressions of L2, sigma and eta into symbols of L2, sigma and eta
    substitutions = [
        (sigma, sigma_sym), # (B.4.3)
        (eta, eta_sym), # (B.4.4)
        (L2, L2_sym), # (B.4.1)
        (sp.sqrt((e2 - 1) * (sigma_sym**2 + eta_sym**2 - 1)), zeta), # derived from (B.4.5)
        ((sigma_sym**2 + eta_sym**2 - 1), -zeta**2 / (1 - e2)), # derived from (B.4.5)
    ]
    g_diff_simplified = g_diff.subs(substitutions).simplify()
    render(g_diff_simplified, name=" ")

    print(f"tan(f_2) * (1 - e^2) = {tan_f2 * (1 - e2_earth):.3e} << 1")
    print("Hence we could assume that $|g'(zeta)| < 1$, hence will converge by fixed-point iteration.")

fixed_point_analysis()

def question_b6():
    N = 4 # Chosen specifically so that $\alpha = \pi/2$ is one of the alphas
    t = -0.87

    alpha, zeta, iters = compute_umbral_boundary(4, -0.87)
    L2 = l2(t) - zeta * tan_f2
    sigma = x(t) + L2 * np.cos(alpha)
    eta = y(t) + L2 * np.sin(alpha)

    print()
    print(f"Solution found in {iters} iterations.")
    print(f"Besselian coordinates of umbral cone with N = {N}:")
    print(f"alpha: \t{alpha}")
    print(f"sigma: \t{sigma}")
    print(f"eta: \t{eta}")
    print(f"zeta: \t{zeta}")
    print()
    print(r"Therefore we found that at $\alpha = \pi/2$ of the Umbral cone, the Besselian coordinates are:")
    print(f"(sigma, eta, zeta) = ({sigma[1]: .5f}, {eta[1]: .5f}, {zeta[1]: .5f})")

question_b6()


# ***
# ### Question B.7.
#
# > Convert the Besselian coordinates computed in Question B.6 for $N = 100$,
# > $t - t_0 = -0.87$ to geographical latitude and longitude using your code
# > from Question B.3. Plot the resulting closed loop on a map of southern
# > England with `cartopy`, along with the centre line. Does the shadow
# > encompass London (geographical latitude $51.5074^{\circ}$, longitude
# > $-0.1278^{\circ}$)?
# >
# > [3 marks]

# References:
# - https://commons.wikimedia.org/wiki/File:Map_of_Southern_England_with_settlements_and_traditional_counties.png to see what constitutes Southern England,
# - https://www.mapsofworld.com/lat_long/united-kingdom-lat-long.html to get the coordinates for the boundary of Southern England,

def question_b7():
    N = 100
    t = -0.87
    # Compute alpha and zeta represending the Umbral boundary points
    # given the parameters provided in the question
    alphas, zetas, iters = compute_umbral_boundary(N, t)
    print(f"Solution found in {iters} iterations.")

    # Calculate the Besselian coordinates of these boundary points
    L2 = l2(t) - zetas * tan_f2
    sigmas = (x(t) + L2 * np.cos(alphas)).reshape(-1, 1)
    etas = (y(t) + L2 * np.sin(alphas)).reshape(-1, 1)
    ts = np.array([t])

    # Convert the Besselian coordinates into latitude and longitude in degrees
    latitudes_umbral, longitudes_umbral = compute_latitude_longitude(sigmas, etas, ts)
    latitudes_umbral = latitudes_umbral.squeeze()
    longitudes_umbral = longitudes_umbral.squeeze()

    # I want matplotlib to close the loop showing the Umbral cone,
    # so append first point to the end of the list as well
    latitudes_umbral = np.append(latitudes_umbral, latitudes_umbral[0])
    longitudes_umbral = np.append(longitudes_umbral, longitudes_umbral[0])

    # Get latitude and longitude coordinates for centre line
    latitudes, longitudes = get_centre_line_coords(start_time, end_time)

    fig = plt.figure(figsize=(12, 7))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Mercator())
    # Using 10m scale here as the map is much more scaled up compared to in Question B5.
    ax.add_feature(cfeature.LAND.with_scale("10m"))
    ax.add_feature(cfeature.OCEAN.with_scale("10m"))
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"))
    # Coordinates obtained from https://www.mapsofworld.com/lat_long/united-kingdom-lat-long.html
    # to cover Southern England, as well as the full Umbral Shadow
    ax.set_extent([-6, 2, 49, 54])

    # Plot the centre line of the eclipse
    ax.plot(longitudes, latitudes, transform=ccrs.Geodetic(), color="red", label="Centre line")

    # Draw the umbral cone
    ax.plot(longitudes_umbral, latitudes_umbral, transform=ccrs.Geodetic(), color="blue")
    # Also shade the area in the umbral cone
    ax.fill(longitudes_umbral, latitudes_umbral, transform=ccrs.Geodetic(), color="blue", alpha=0.3, label="Umbral shadow")

    # Mark the position of London in the map
    # Using linestyle None simply for aesthetic reason, so that no line
    # appears in the legend.
    london_coords = (-0.1278, 51.5074) # as provided
    ax.plot(london_coords[0], london_coords[1], transform=ccrs.Geodetic(), color="green", marker="o", linestyle="None", label="London")

    # Compute the time which we drew the umbral shadow
    shadow_time = t0 + TimeDelta(t * 3600, format="sec")
    
    # Add chart elements
    ax.gridlines(draw_labels=True, dms=False, x_inline=False, y_inline=False)
    plt.legend()
    plt.title(f"1715 Solar Eclipse Centre Line and Umbral Shadow\nover Southern England at {shadow_time.strftime("%H:%M:%S")} TDT")

question_b7()

# We can observe visually that London does indeed lie within the Umbral cone.

# We can check that this map looks fairly accurate comparing to the map made by Edmond Halley
# https://commons.wikimedia.org/wiki/File:Solar_eclipse_1715May03_Halley_map.png
