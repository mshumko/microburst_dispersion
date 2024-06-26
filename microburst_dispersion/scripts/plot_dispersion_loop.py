from datetime import datetime, date
import dateutil.parser
import string

import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.dates
from matplotlib.ticker import FuncFormatter
import matplotlib.gridspec as gridspec
import pandas as pd
import numpy as np

import microburst_dispersion
import microburst_dispersion.firebird
import microburst_dispersion.fit_duration


class Dispersion_Summary():
    """
    Make a series of summary plots for the microbursts.
    """
    def __init__(
            self, sc_id:int, catalog_version:int=5, plot_window_s:float=1, fit_interval_s:float=0.3, annotate_fit:bool=False
            ) -> None:
        self.sc_id = sc_id
        self.catalog_version = catalog_version
        self.fit_interval_s = pd.Timedelta(seconds=fit_interval_s)
        self.plot_window_s = pd.Timedelta(seconds=plot_window_s)
        self.current_date = date.min
        self.channels:list=np.arange(6)
        self._plot_colors = np.array(['k', 'k', 'k', 'k', 'k', 'k'])
        self.full_ylabels = True
        self.annotate_fit = annotate_fit

        catalog_name = f'FU{sc_id}_microburst_catalog_{str(catalog_version).zfill(2)}.csv'
        catalog_path = microburst_dispersion.config['here'].parent / 'data' / catalog_name
        self.catalog = pd.read_csv(catalog_path)
        self.catalog['Time'] = pd.to_datetime(self.catalog['Time'])
        return

    def loop(self, debug=True):
        self.t0_keys = [f't0_{channel}' for channel in self.channels]
        if not debug:
            save_dir = microburst_dispersion.config['here'].parent / 'plots' / 'dispersion_summary'
            if not save_dir.exists():
                save_dir.mkdir()
                print(f'Made {save_dir} directory')

        for i, row in self.catalog.iterrows():
            print(f'Progress: {i}/{self.catalog.shape[0]}')
            if np.any(pd.isnull(row[self.t0_keys].to_numpy())):
                continue  # The fit failed somewhere.
            self.microburst_info = row
            current_date = date.min
            if current_date != row['Time'].date():
                self.hr = microburst_dispersion.firebird.Hires(self.sc_id, row['Time']).load()
                self.cadence_s = float(self.hr.attrs["CADENCE"])
                self.cadence_ms = 1000*self.cadence_s
                self.center_energy, self.energy_range, self.energy_range_array = self.get_energy_channels()
                current_date = row['Time'].date()
            time_range = (
                row['Time']-self.plot_window_s/2, 
                row['Time']+self.plot_window_s/2
                )
            self.plot_idt = np.where(
                (self.hr['Time'] > time_range[0]) & (self.hr['Time'] < time_range[1])
                )[0]
            self._create_subplots()
            self._plot_hr()
            self._plot_fit()
            self._annotate_location()
            self._plot_dispersion(self.ax[-1])
            self.ax[0].set_title(f'FU{self.sc_id} Microburst Dispersion\n{self.microburst_info["Time"]:%F %T}')
            self._format_times(self.ax[-2])
            self.ax[-2].set_xlabel('Time [HH:MM:SS]')

            # plt.tight_layout()
            if debug:
                plt.show()
            else: 
                save_name = f'{self.microburst_info["Time"]:%Y%m%d_%H%M%S}_fu{self.sc_id}_microburst_dispersion.png'
                plt.savefig(save_dir/save_name)
                plt.close()
        return
    
    def get_energy_channels(self):
        center_energy = np.array([float(s.split()[0].replace('>', '')) 
                              for s in np.array(self.hr.attrs['Col_counts']['ELEMENT_LABELS'])])
        center_energy = center_energy[self.channels]
        energy_range = np.array(self.hr.attrs['Col_counts']['ENERGY_RANGES'])
        energy_range = energy_range[self.channels]

        # Turn the energy range from strings to a (n_channels, 2) array of energy channel bounds.
        energy_range_array = np.zeros((len(self.channels), 2), dtype=float)
        for i in range(len(self.channels)):
            _energy_range_row = energy_range[i].replace(' ', '')
            if '>' in _energy_range_row:
                # Special case with integral channel.
                _lower_intergral_energy= int(float(_energy_range_row[:-3].replace(">", "")))
                # If memory serves me right, 2 MeV is effective upper bound on the integral energy channel
                energy_range_array[i, :] = [_lower_intergral_energy, 2000]
            else:
                # DIfferential channels.
                energy_range_array[i, :] = np.array(_energy_range_row[:-3].split('-'), dtype=float).astype(int)
        return center_energy, energy_range, energy_range_array.astype(int)
    
    def _create_subplots(self):
        """"
        Create empty subplots for the HiRes line plots and dispersion scatter plot.
        """
        # I want to adjust the hspace for the HiRes line subplots and the dispersion 
        # subplot separately so I created multiple nested gridspecs.
        # See https://stackoverflow.com/a/31485288 for inspiration
        outer_gridspec = gridspec.GridSpec(2, 1, height_ratios=[len(self.channels), 1], 
                                           top=0.94, left=0.16, right=0.958, bottom=0.055, hspace=0.15) 
        inner_gs1 = gridspec.GridSpecFromSubplotSpec(len(self.channels), 1, subplot_spec=outer_gridspec[0], hspace=0.05)
        inner_gs2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_gridspec[1])

        self.fig = plt.figure(figsize=(6, 8))
        self.ax = [None]*(len(self.channels)+1)
        for i in range(len(self.channels)):
            if i == 0:
                self.ax[i] = self.fig.add_subplot(inner_gs1[i, 0])
            else:
                self.ax[i] = self.fig.add_subplot(inner_gs1[i, 0], sharex=self.ax[0])
            if i < len(self.channels)-1:
                self.ax[i].get_xaxis().set_visible(False)
        self.ax[-1] = self.fig.add_subplot(inner_gs2[0, 0])
        for i, (ax_i, color) in enumerate(zip(self.ax, self._plot_colors)):
            ax_i.text(0, 0.99, f'({string.ascii_uppercase[i]})', va='top', 
                      transform=ax_i.transAxes, weight='normal', fontsize=13, color=color)

        x_offset = -0.17
        # -1 since it is in reference to ax[0] and 0.5 to put it in the middle
        text_y_center = -((len(self.channels)-1)/2)+0.5
        plt.annotate("Energy",
               xy=(x_offset, 1), xycoords=self.ax[0].transAxes,
               xytext=(x_offset, text_y_center), textcoords=self.ax[0].transAxes,
               arrowprops=dict(arrowstyle="-|>", lw=1, color='black'), rotation='vertical', 
               ha='center', va='center', fontsize=15)
        plt.annotate("Energy",
               xy=(x_offset, 0), xycoords=self.ax[-2].transAxes,
               xytext=(x_offset, text_y_center), textcoords=self.ax[0].transAxes,
               arrowprops=dict(arrowstyle="-", lw=1), rotation='vertical', 
               ha='center', va='center', fontsize=15)
        return
    
    def _format_times(self, ax):
        locator=matplotlib.ticker.MaxNLocator(nbins=5)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(FuncFormatter(self.format_fn))
        return

    def format_fn(self, tick_val, _):
        tick_time = matplotlib.dates.num2date(tick_val).replace(tzinfo=None)
        return tick_time.strftime('%T.%f')[:-5]
    
    def _plot_hr(self):
        """
        Plot the Hires channels in each subplot.
        """
        # Plot HiRes
        for i, (ax_i, color, channel) in enumerate(zip(self.ax[:-1][::-1], self._plot_colors, self.channels)):
            ax_i.step(
                self.hr['Time'][self.plot_idt], 
                self.hr['Col_counts'][self.plot_idt, channel]/self.cadence_s, c=color, where='mid'
                )
            
            if self.full_ylabels:
                ax_i.set_ylabel(f'[counts/s]')
                ax_i.text(0.07, 0.99, f'{self.energy_range_array[i, 0]} - {self.energy_range_array[i, 1]} keV', va='top', 
                      transform=ax_i.transAxes, weight='normal', color=color, fontsize=13)
            else:
                ax_i.set_ylabel(f'{channel=}')
            max_counts = np.max(self.hr['Col_counts'][self.plot_idt, channel])/self.cadence_s
            ax_i.set_ylim(0, 1.2*max_counts)
        return
    
    def _annotate_location(self):
        lat_str = f'${{{round(self.microburst_info["Lat"])}}}^{{\circ}}$'
        lon_str = f'${round(self.microburst_info["Lon"])}^{{\circ}}$'
        s = (
            f'L={round(self.microburst_info["McIlwainL"], 1)}\n'
            f'MLT={round(self.microburst_info["MLT"], 1)}\n'
            f'(lat,lon)=({lat_str},{lon_str})'
            )
        self.ax[0].text(0.67, 1, s, va='top', transform=self.ax[0].transAxes, color='k')
        return
    
    def _plot_fit(self):
        """
        Plot the Gaussian fit + linear trend
        """
        # Plot Fit
        data_time_array = self.hr['Time'][self.plot_idt]
        time_array = pd.date_range(start=data_time_array[0], end=data_time_array[-1], periods=1000)
        current_date = time_array[0].floor('d')
        x_data_seconds = (time_array-current_date).total_seconds()

        for ax_i, color, channel in zip(self.ax[:-1][::-1], self._plot_colors, self.channels):
            fit_bounds = (
                self.microburst_info['Time']-self.fit_interval_s/2,
                self.microburst_info['Time']+self.fit_interval_s/2
            )
            ax_i.axvspan(*fit_bounds, color='grey', alpha=0.5)
            ax_i.axvline(self.microburst_info['Time'], color='k', ls=':')
            
            popt = np.nan*np.zeros(5)
            popt[1] = (dateutil.parser.parse(self.microburst_info[f't0_{channel}']) - current_date).total_seconds()
            popt[2] = self.microburst_info[f'fwhm_{channel}']/2.355 # Convert the Gaussian FWHM to std
            popt[0] = self.microburst_info[f'A_{channel}']
            popt[3] = self.microburst_info[f'y_int_{channel}']
            popt[4] = self.microburst_info[f'slope_{channel}']
            if np.isnan(popt[0]):  # Plot just the data if the fit failed.
                continue

            gaus_y = microburst_dispersion.fit_duration.Duration.gaus_lin_function(x_data_seconds, *popt)
            ax_i.plot(time_array, gaus_y/self.cadence_s, c=color, ls='--')

            if self.annotate_fit:
                fit_params=(
                    f"FWHM={round(self.microburst_info[f'fwhm_{channel}'], 2)} [s]\n"
                    f"R^2 = {round(self.microburst_info[f'r2_{channel}'], 2)}\n"
                    f"adj_R^2 = {round(self.microburst_info[f'adj_r2_{channel}'], 2)}\n"
                )
                ax_i.text(0.01, 0.87, fit_params, va='top', transform=ax_i.transAxes, color=color)
        return
    
    def _plot_dispersion(self, ax):
        self._get_dispersion()
        ax.errorbar(self.center_energy, self.t0_diff_ms, c='k', marker='.', 
            yerr=self.yerrs, xerr=self.xerrs, capsize=2, ls='None')
        # max_abs_lim = 1.1*np.max(np.abs(ax.get_ylim()))
        ax.set_ylim(-50, 50)
        ax.axhline(c='k', ls='--')
        ax.set(xlabel='Energy [keV]', ylabel='Peak time delay [ms]\n(ch[N]-ch[0])')

        locator=matplotlib.ticker.FixedLocator(np.linspace(-50, 50, num=5))
        ax.yaxis.set_major_locator(locator)
        return
    
    def _get_dispersion(self):
        # Time differences with respect to channel 0
        t0 = [dateutil.parser.parse(self.microburst_info[t0_key]) for t0_key in self.t0_keys]
        self.t0_diff_ms = [1E3*(t0_i - t0[0]).total_seconds() for t0_i in t0]

        self.xerrs = [xerr for xerr in self.energy_range]
        self.xerrs = [xerr.replace('keV', '').replace('>', '').split('-') for xerr in self.xerrs]

        if 5 in self.channels:  # Integral channel special case
            self.xerrs[-1] = [None, None]
        self.xerrs = np.array(self.xerrs).astype(float).T - self.center_energy
        self.xerrs = np.abs(self.xerrs)
        self.yerrs = self.cadence_ms/2
        return

if __name__ == '__main__':
    catalog_version:int=5
    plot_window_s = 1

    for sc_id in [4]:
        d = Dispersion_Summary(
            sc_id, 
            catalog_version=catalog_version, 
            plot_window_s=plot_window_s,
            annotate_fit=False
            )
        d.loop(debug=False)