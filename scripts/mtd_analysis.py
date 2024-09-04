import pandas as pd
import numpy as np
from lmfit import Model


def integrate(xvals, yvals, xlower, xupper):
    """
    Integrate y over x from xlower to xupper.
    Use trapz to integrate over points closest to xlower, xupper.
    the +1 to idx_max is for numpy half-open indexing.
    """

    idx_min = np.argmin(np.abs(xvals - xlower))
    idx_max = np.argmin(np.abs(xvals - xupper)) + 1
    result = np.trapz(yvals[idx_min:idx_max], x=xvals[idx_min:idx_max])
    return result


def mtd_line_binding(
    fes_fname: str,
    spring_constant: float,
    temperature: float = 300,
    radius: float = 0.0,
) -> float:

    # Standard volume
    STANDARD_V0 = 1661 / 1000  # nm^-3
    MOLAR_GAS_CONSTANT_R = 8.31446261815324  # J / mol / K

    kBT = MOLAR_GAS_CONSTANT_R / 1000 * temperature  # kJ / mol

    # integral limits
    b0, b1, b2 = -0.5, 0.3, 1.4  # nm
    unbound_length = b2 - b1  # nm

    print("---------------------- Binding Energy ---------------------")
    print(f"FES Filename   : {fes_fname:<20}")
    print(f"Spring constant: {spring_constant:10.10f} kJ / mol / nm^2")
    print(f"Cylinder radius: {radius:10.10f} nm")
    print(f"Temperature    : {temperature:10.10f} K")
    print(f"  kBT          : {kBT:10.10f} kJ/mol")
    print("Integral limits:")
    print(f" Bound         : {b0:10.10f} - {b1:10.10f} nm")
    print(f" Unbound       : {b1:10.10f} - {b2:10.10f} nm")
    print(f"Standard volume: {STANDARD_V0:10.10f} nm^3")

    # 1D free energy on a line
    data = pd.read_csv(fes_fname, header=None, sep="\\s+", comment="#")
    xval = data.iloc[:, 0]
    yval = data.iloc[:, 1]

    mean_unbound = np.mean(yval.loc[(xval > 0.8) & (xval < b2)])

    write = pd.DataFrame({"z": xval, "dG": yval - mean_unbound})
    with open("aligned_fes.dat", "w") as f:
        f.write("# z dG")
    write.to_csv("aligned_fes.dat", index=False, header=False, mode="a", sep=" ")
    print("Potential energy surface aligned to unbound region in file aligned_fes.dat")

    yval = np.exp(-yval / kBT)

    unbound_volume = unbound_length * (
        2
        * np.pi
        * (
            0.5 * radius**2
            + radius * np.sqrt(np.pi * kBT / 2 / spring_constant)
            + kBT / spring_constant
        )
    )
    print(f"Unbound volume : {unbound_volume:10.10f} nm^3")

    bound_integral = integrate(xval, yval, b0, b1)
    unbound_integral = integrate(xval, yval, b1, b2)

    binding_dG = -kBT * np.log(unbound_integral / bound_integral)
    volume_dG = kBT * np.log(unbound_volume / STANDARD_V0)
    total_dG = volume_dG + binding_dG
    print("-----------------------------------------------------------")
    print(f"Binding free energy (dG)        : {binding_dG:10.2f} kJ/mol")
    print(f"Volume correction   (dG)        : {volume_dG:10.2f} kJ/mol")
    print(f"Total corrected free energy (dG): {total_dG:10.2f} kJ/mol")
    print("-----------------------------------------------------------")
    return total_dG


def estimate_bound_volume_correction(
    spring_constant: float, temperature: float = 300, radius: float = 0.0
):
    print("--------------------- Bound correction --------------------")

    def global_pes(x):
        return -30 * np.exp(-((x) ** 2 / (2 * (0.1**2)))) + 15 * np.exp(
            -((x - 0.3) ** 2 / (2 * (0.1**2)))
        )

    MOLAR_GAS_CONSTANT_R = 8.31446261815324  # J / mol / K
    kBT = MOLAR_GAS_CONSTANT_R / 1000 * temperature  # kJ / mol

    x = np.linspace(-1, 2, num=10000)
    y = global_pes(x)
    offset = np.min(y)
    y -= offset  # translate minimum to zero

    b0, b1 = -0.15, 0.15  # nm
    print("Fit limits:")
    print(f" Lower         : {b0:10.10f} nm")
    print(f" Upper         : {b1:10.10f} nm")
    indc = np.argwhere((x > b0) & (x < b1)).flatten()
    x = x[indc]
    y = y[indc]

    def harmonic(x, k):
        return 0.5 * k * x**2

    harmonic_model = Model(harmonic)
    params = harmonic_model.make_params(k=1000)
    fit = harmonic_model.fit(data=y, params=params, x=x)
    fit_spring_constant = fit.params["k"].value
    print(f"Fitted spring constant: {fit_spring_constant:10.10f} kJ / mol / nm^2")
    fit_y = harmonic(x, k=fit_spring_constant)
    write = pd.DataFrame({"z": x, "dG": fit_y + offset})
    with open("fit_harmonic.dat", "w") as f:
        f.write("# z dG")
    write.to_csv("fit_harmonic.dat", index=False, header=False, mode="a", sep=" ")
    print("Harmonic fit written to fit_harmonic.dat")

    v0 = (
        4
        * np.pi
        * (kBT / fit_spring_constant)
        * np.sqrt((np.pi * kBT) / (2 * fit_spring_constant))
    )
    print(f"Estimated unrestrained bound volume: {v0:10.10f} nm^3")
    restrained_volume = (
        (
            2
            * np.pi
            * (
                0.5 * radius**2
                + radius * np.sqrt((np.pi * kBT) / (2 * spring_constant))
                + kBT / spring_constant
            )
        )
        * (np.sqrt((np.pi * kBT) / (2 * fit_spring_constant)))
        * 2
    )
    volume_dG = -kBT * np.log(restrained_volume / v0)
    if volume_dG < 0.0:
        volume_dG = 0.0
    print(f"Estimated restrained bound volume: {restrained_volume:10.10f} nm^3")
    print("-----------------------------------------------------------")
    print(f"Estimated volume correction free energy (dG): {volume_dG:10.2f} kJ/mol")
    print("-----------------------------------------------------------")
    return volume_dG


if __name__ == "__main__":
    import sys

    fes_fname = str(sys.argv[1])
    spring_constant = float(sys.argv[2])  # kJ / mol / nm^2
    temperature = float(sys.argv[3])  # Kelvin
    radius = float(sys.argv[4])  # nm
    be = mtd_line_binding(
        fes_fname=fes_fname,
        spring_constant=spring_constant,
        temperature=temperature,
        radius=radius,
    )
    vc = estimate_bound_volume_correction(
        spring_constant=spring_constant, temperature=temperature, radius=radius
    )
    print(f"Final corrected binding free energy: {be + vc: 10.2f} kJ/mol")
