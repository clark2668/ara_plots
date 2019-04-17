from limits import LimitFigure, ara_energies, ara_available, ara_projected, efficiencies


if __name__=="__main__":
    figure = LimitFigure(e_power=1, xlims=(1e6, 1e11), ylims=(1e-19, 2e-14), font_size=16, tick_size=14)
    figure.build_base_plot('ara')

    ara_available=ara_available*efficiencies
    ara_projected=ara_projected*efficiencies

    energy_a21yr, flux_a21yr, _, _ = figure.get_data('sensitivities/ara_ming.txt')
    energy_a24yr, flux_a214r, _, _ = figure.get_data('sensitivities/ara_projected_eff.txt')
    _plt_a21yr, = figure.ax.plot(energy_a21yr, flux_a21yr, color='#2288AA',linestyle=None,label='ARA2 Measured (2013)',linewidth=3)
    _plt_a24yr, = figure.ax.plot(energy_a24yr, flux_a214r, color='#2288AA',linestyle='--',label='ARA2 Expected (2013-2016)',linewidth=3)
 
    figure.custom_limits.append(_plt_a21yr)
    figure.custom_limits.append(_plt_a24yr)

    figure.add_limit('ARA', ara_energies, ara_available,
                    stations=2.44, years=1, color='black', linestyle=None,
                    label='ARA5 Projected (2012-2019)')
    figure.add_limit('ARA', ara_energies, ara_projected,
                    stations=2.44, years=1, color='black', linestyle='--',
                    label='ARA5 Projected (2012-2022)')


    figure.show(legend_size=10, save_name='e1_sensitivity_plot_ara.pdf')


    # figure = LimitFigure(e_power=1, xlims=(1e6, 1e11), ylims=(1e-19, 2e-14), font_size=16, tick_size=14)
    # figure.build_base_plot('ara_src')
    # figure.add_limit('ARA', ara_energies, ara_available,
    #                 stations=2.44, years=1, color='black', linestyle=None,
    #                 label='Available Data (2012-2019)')
    # figure.add_limit('ARA', ara_energies, ara_projected,
    #                 stations=2.44, years=1, color='black', linestyle='--',
    #                 label='Projected Data (2012-2022)')
    # # figure.title("Trigger Level Sensitivities")
    # figure.show(legend_size=10, save_name='e1_source_plot_ara.pdf')