import math
import numpy as np
import matplotlib.pyplot as plt


class BGC_MultipleDosing:
    def __init__(self, properties):
        self.props = properties
        self.dosing_profile = None  # [{dose:<amt>, tau=<amt>, tobs=[], conc=[]}, {profile_2}]
        self.t_obs = None

    def set_observation_times(self, min, max, interval):
        # set observation times for class
        self.t_obs = list(np.arange(min, max, interval))

    # concentration time values
    def Ct_MultDose(self, D0, t, tau):
        n_doses = np.ceil(t / tau)
        tsld = t - (n_doses - 1) * tau

        first_term = (self.props['F'] * D0 * self.props['ka']) / (
                    self.props['Vd'] * (self.props['ka'] - self.props['kel']))

        f_exp_kel = (1 - math.exp(-self.props['kel'] * tau * n_doses)) / (1 - math.exp(-self.props['kel'] * tau))
        first_exp = f_exp_kel * math.exp(-self.props['kel'] * tsld)

        f_exp_ka = (1 - math.exp(-self.props['ka'] * tau * n_doses)) / (1 - math.exp(-self.props['ka'] * tau))
        second_exp = f_exp_ka * math.exp(-self.props['ka'] * tsld)

        # print("t={0}, dose={1}, first_term={2}, f_exp_kel={3}, first_exp={4}, f_exp_ka={5}, second_exp={6}".
        #       format(t, n_doses, first_term, f_exp_kel, first_exp, f_exp_ka, second_exp))

        return first_term * (first_exp - second_exp)

    def simulate_profile(self, dose, tau, ret=False):
        # simulate concentration time profile
        if self.t_obs is None:
            print("Please set observation times before simulating profile.")
            raise Exception

        c_profile = []
        for t in self.t_obs:
            c_t = self.Ct_MultDose(D0=dose, t=t, tau=tau)
            c_profile.append(c_t)

        d_prof_dict = dict(dose=dose, tau=tau, tobs=self.t_obs, conc=c_profile)
        d_prof_dict = self.percent_in_efficacy_window(d_prof_dict)
        # print(d_prof_dict)
        if ret:
            return [d_prof_dict]
        else:
            self.dosing_profile = [d_prof_dict]

    def add_profile(self, dose, tau):
        # add additional simulation profile or create profile
        if self.dosing_profile is None:
            self.simulate_profile(dose=dose, tau=tau)
        else:
            temp_profile = self.simulate_profile(dose=dose, tau=tau, ret=True)
            # print("temp_profile={0}".format(temp_profile))
            self.dosing_profile.extend(temp_profile)

    def percent_in_efficacy_window(self, profile):
        # calculate the percent of time points in the simulated window
        efficacy_in = [i for i in profile['conc'] if ((i <= self.props['max_eff']) &
                                                      (i >= self.props['min_eff']))]
        # add percent in efficacy window to profile
        profile['f_win'] = round(len(efficacy_in) / len(profile['conc']), 2)

        return profile

    def plot_class_profiles(self):
        # plot concentration-time profiles for all entered treatment regimens
        if self.dosing_profile is None:
            print("Cannot plot empty dosing profile.")
            raise Exception
        # iterate dosing profiles
        fig, ax = plt.subplots()
        ax.axhline(y=self.props['min_eff'], color="black", linestyle="--", linewidth=0.75)
        ax.axhline(y=self.props['max_eff'], color="black", linestyle="--", linewidth=0.75)
        ax.fill_between(self.t_obs, self.props['min_eff'], self.props['max_eff'], color='gray', alpha=0.2)
        ax.set_xlabel("Time (hr)")
        ax.set_ylabel("Sample Concentration")
        for p in self.dosing_profile:
            p_win = round(p['f_win'] * 100, 1)
            ax.plot(p['tobs'], p['conc'], label="dose={0}, tau={1}, p_win={2}%".format(p['dose'], p['tau'], p_win))
        ax.legend(loc="center right", bbox_to_anchor=(1.7, 0.5), ncol=1, fancybox=True,
                  shadow=True)
        fig.show()
