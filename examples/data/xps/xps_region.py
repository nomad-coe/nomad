"""Fitting functions for XPS spectra"""
from typing import Optional
from dataclasses import dataclass
import numpy as np
from numpy.linalg import norm
import plotly.graph_objects as go
import pandas as pd
from h5py import File as H5File
from lmfit import Model, CompositeModel


@dataclass
class XPSRegion:
    """An XPS region representation"""

    binding_energy: np.ndarray
    counts: np.ndarray
    counts_err: np.ndarray
    baseline: Optional[np.ndarray] = None
    _fit_region: slice = slice(None, None)
    _fit_mod: Optional[Model] = None
    fit_result: Optional[Model] = None

    @staticmethod
    def load(filename: str, entry: str) -> "XPSRegion":
        """Load from a NeXus file.

        Args:
            filename (str): The NeXus file name to load.
            entry (str):
                The entry from which to load data.
                Should be the name of an NXentry within the NeXus file.

        Returns:
            XPSRegion: The XPSRegion class with the loaded data.
        """
        with H5File(filename, "r") as xps_file:
            binding_energy = xps_file[f"/{entry}/data/energy"][:]
            cps = xps_file[f"/{entry}/data/data"][:]
            cps_err = xps_file[f"/{entry}/data/data_errors"][:]

        return XPSRegion(binding_energy=binding_energy, counts=cps, counts_err=cps_err)

    def fit_region(self, start: int, stop: int) -> "XPSRegion":
        """Select a fit region within this XPSregion by x-axis value.
        The fit region is always selected between start and stop, regardless of their order.
        Both points are included in the region,
        hence the actual selected region may be a little larger.

        Args:
            start (int): The start ot the region.
            stop (int): The end of the region.

        Returns:
            XPSRegion: This class
        """
        region = np.argwhere(
            (self.binding_energy >= start) & (self.binding_energy <= stop)
        )

        self._fit_region = slice(region[0, 0], region[-1, 0], 1)
        return self

    def fit_model(self, model: Model) -> "XPSRegion":
        """Supply a fit model to fit this xps region.

        Args:
            model (lmfit.Model): The lmfit model to use.

        Returns:
            XPSRegion: This class
        """
        self._fit_mod = model
        return self

    def fit(self, *args, **kwargs) -> "XPSRegion":
        """Perform a fit of the data. You need to define a fit_model first and
        execute a baseline correction before using this method.

        Raises:
            ValueError: If no fit model is provided or the baseline has not been corrected.

        Returns:
            XPSRegion: This class
        """
        if self._fit_mod is None:
            raise ValueError("You need to provide a fit model before performing a fit.")

        if self.baseline is None:
            raise ValueError(
                "You need to perform a baseline correction before using this method."
            )

        self.fit_result = self._fit_mod.fit(
            (
                self.counts[self._fit_region] - self.baseline
                if self._fit_region
                else self.counts[self._fit_region]
            ),
            *args,
            x=self.binding_energy[self._fit_region],
            weights=1 / self.counts_err[self._fit_region],
            **kwargs,
        )

        return self

    def calc_baseline(self, bg_type: str = "shirley") -> "XPSRegion":
        """Calculate the baseline for this xps spectrum in the given region.

        Args:
            bg_type (str, optional): The background type. Defaults to "shirley".

        Raises:
            ValueError: If the bg_type is unsupported.

        Returns:
            XPSRegion: This class
        """
        baselines = {"shirley": shirley_baseline}
        if bg_type not in baselines:
            raise ValueError(f"Unsupported baseline type {bg_type}.")

        self.baseline = baselines[bg_type](
            self.binding_energy[self._fit_region], self.counts[self._fit_region]
        )

        return self

    def peak_property(self, prop: str) -> pd.DataFrame:
        """Generates a dataframe with values for a property `prop` of a fitting model.

        Args:
            prop (str):
                The name of the property to deduce for the peaks,
                e.g. `center` for the peak center for gaussian or lorentzian shapes.

        Raises:
            ValueError:  Thrown if no prior fit is performed.

        Returns:
            pd.DataFrame: A pandas DataFrame containing the peak property for each peak.
        """
        if not self.fit_result:
            raise ValueError("You need to perform a fit first.")

        props = pd.DataFrame()
        for prefix in map(lambda x: x.prefix, self.fit_result.components):
            if f"{prefix}{prop}" not in self.fit_result.params:
                continue

            props = pd.concat(
                [
                    props,
                    pd.DataFrame(
                        {prop: self.fit_result.params.get(f"{prefix}{prop}").value},
                        index=[prefix.rstrip("_")],
                    ),
                ]
            )

        return props

    def peak_areas(self, region_only=False) -> pd.DataFrame:
        """Calculates the peak areas of the given fit models peaks.

        Args:
            region_only (bool, optional):
                Set true if only the area inside the set region should be consider.
                Defaults to False.

        Raises:
            ValueError: Thrown if no prior fit is performed.

        Returns:
            pandas.DataFrame: A pandas DataFrame containing the peak areas.
        """
        if not self.fit_result:
            raise ValueError("You need to perform a fit first.")

        areas = pd.DataFrame()
        if region_only:
            peaks = self.fit_result.eval_components(
                x=self.binding_energy[self._fit_region]
            )
        else:
            peaks = self.fit_result.eval_components(x=self.binding_energy)
        for prefix in peaks:
            areas = pd.concat(
                [
                    areas,
                    pd.DataFrame(
                        {"Area": sum(peaks[prefix])},
                        index=[prefix.rstrip("_")],
                    ),
                ]
            )
        return areas

    def plot_residual(self):
        """Plot the fit residual"""
        if not self.fit_result:
            raise ValueError("You need to perform a fit first.")

        fig = go.Figure(
            data=go.Scatter(
                name="Residual",
                x=self.binding_energy[self._fit_region],
                y=self.fit_result.residual,
            )
        )
        fig.update_xaxes(title="Binding energy (eV)")
        fig.update_yaxes(title="Residual")
        fig.layout.xaxis.autorange = "reversed"

        fig.show()

    def plot(self):
        """Plot the xps region"""
        fig = go.Figure(
            data=go.Scatter(
                name="Measurement",
                x=self.binding_energy,
                y=self.counts,
                error_y=dict(
                    type="data",  # value of error bar given in data coordinates
                    array=self.counts_err,
                    visible=True,
                ),
            )
        )
        if self._fit_region.start is not None:
            fig.add_vline(
                self.binding_energy[self._fit_region.start], line_color="lightgrey"
            )
        if self._fit_region.stop is not None:
            fig.add_vline(
                self.binding_energy[self._fit_region.stop], line_color="lightgrey"
            )
        if self.baseline is not None:
            fig.add_trace(
                go.Scatter(
                    name="Baseline",
                    x=self.binding_energy[self._fit_region],
                    y=self.baseline,
                ),
            )
        if self.fit_result is not None:
            fig.add_trace(
                go.Scatter(
                    name="Fit",
                    x=self.binding_energy[self._fit_region],
                    y=self.fit_result.best_fit + self.baseline,
                ),
            )
        if isinstance(self._fit_mod, CompositeModel):
            peaks = self.fit_result.eval_components(
                x=self.binding_energy[self._fit_region]
            )
            for prefix in peaks:
                fig.add_trace(
                    go.Scatter(
                        name=prefix.rstrip("_"),
                        x=self.binding_energy[self._fit_region],
                        y=peaks[prefix] + self.baseline,
                    )
                )

        fig.update_xaxes(title="Binding energy (eV)")
        fig.update_yaxes(title="CPS")
        fig.layout.xaxis.autorange = "reversed"

        fig.show()


# pylint: disable=invalid-name
def shirley_baseline(
    x: np.ndarray, y: np.ndarray, tol: float = 1e-5, maxit: int = 10
) -> np.ndarray:
    """Calculate the shirley background according to the Sherwood method.

    Args:
        x (np.ndarray): The x-axis on which to calculate the shirley baseline.
        y (np.ndarray): The y-axis on which to calculate the shirley baseline.
        tol (float, optional):
            The convergence tolerance at which to stop the iteration.
            Defaults to 1e-5.
        maxit (int, optional):
            The maximum iteration after which to stop the iteration.
            Defaults to 10.

    Raises:
        ValueError:
            Is thrown when the arrays have different dimensions,
            are not numpy arrays or are empty.
            Is also thrown when the fit does not converge.

    Returns:
        np.ndarray: The shirley baseline for the x, y dataset.
    """

    if not isinstance(x, np.ndarray) or not isinstance(y, np.ndarray):
        raise ValueError(
            f"Parameters x and y must be of type numpy array, not {type(x)} and {type(y)}"
        )

    if len(x) != len(y):
        raise ValueError("x and y arrays have different dimensions.")

    if not x.any():
        raise ValueError("x-array is empty.")

    if len(x.shape) > 1:
        raise ValueError(
            f"Data arrays must be one-dimensional. Found dimension {x.shape}."
        )

    is_reversed = False
    if x[0] < x[-1]:
        is_reversed = True
        x = x[::-1]
        y = y[::-1]

    background = np.zeros(x.shape)
    background_next = background.copy()

    iters = 0
    while True:
        k = (y[0] - y[-1]) / np.trapz(y - background, x=x)

        for energy in range(len(x)):
            background_next[energy] = k * np.trapz(
                y[energy:] - background[energy:], x=x[energy:]
            )

        diff = norm(background_next - background)
        background = background_next.copy()
        if diff < tol:
            break

        iters += 1
        if iters == maxit:
            raise ValueError(
                "Maximum number of iterations exceeded before convergence."
            )

    if is_reversed:
        return (y[-1] + background)[::-1]
    return y[-1] + background
