# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#                                    Author : Ronan Arraes Jardim Chagas <ronisbr@gmail.com>
#                                      Date : 2023-04-25T13:00:00 UTC-3
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Create a test set for the J2 osculating orbit propagator using a numerical algorith.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using DifferentialEquations
using PrettyTables
using SatelliteToolbox
using StaticArrays

"""
    create_j2_osc_validation_tests()

Create the validation tests for the J2 osculating orbit propagator.
"""
function create_j2_osc_validation_tests()
    # Load gravity model
    # ======================================================================================

    # We are using the constants in the EGM-2008 model
    icgem_egm08 = parse_icgem("../EGM2008.gfc")
    egm08 = create_gravity_model_coefs(icgem_egm08)

    # Inputs
    # ======================================================================================

    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
    orb = KeplerianElements(
        jd₀,
        7190.982e3,
           0.001111,
          98.405 |> deg2rad,
          90     |> deg2rad,
         200     |> deg2rad,
          45     |> deg2rad
    )
    orbp = init_orbit_propagator(Val(:J2osc), orb)

    # Initial osculating value
    # ======================================================================================

    r₀, v₀ = propagate!(orbp, 0)

    # Build the ODE problem
    # ======================================================================================

    p = (; jd₀, egm08)
    u₀ = vcat(r₀, v₀)
    tspan = (0.0, 6000.0)
    prob = ODEProblem(_dynamics, u₀, tspan, p)

    # Solve and return the data
    # ======================================================================================

    sol = solve(
        prob,
        AutoVern7(Rodas5());
        abstol = 1e-14,
        reltol = 1e-14,
        saveat = tspan[1]:600:tspan[2]
    )

    t  = sol.t
    rx = map(x -> x[1], sol.u)
    ry = map(x -> x[2], sol.u)
    rz = map(x -> x[3], sol.u)
    vx = map(x -> x[4], sol.u)
    vy = map(x -> x[5], sol.u)
    vz = map(x -> x[6], sol.u)

    pretty_table(
        hcat(t, rx ./ 1000, ry ./ 1000, rz ./ 1000, vx ./ 1000, vy ./ 1000, vz ./ 1000);
        formatters = ft_printf("%.3f", 2:7),
        header = (
            [
                "Time [s]",
                "Pos. X (TOD)",
                "Pos. Y (TOD)",
                "Pos. Z (TOD)",
                "Vel. X (TOD)",
                "Vel. Y (TOD)",
                "Vel. Z (TOD)"
            ],
            ["", "km", "km", "km", "km / s", "km / s", "km / s"]
        ),
        hlines = [:header],
        vlines = [1]
    )

    return nothing
end

############################################################################################
#                                    Private Functions
############################################################################################

function _dynamics(u, p, t)
    # Unpack the parameters.
    jd₀   = p.jd₀
    egm08 = p.egm08

    # Unpack the state.
    r_tod = SVector(u[1], u[2], u[3])
    v_tod = SVector(u[4], u[5], u[6])

    # Compute the matrix to convert the TOD reference frame to the PEF reference frame.
    D_pef_tod = r_eci_to_ecef(TOD(), PEF(), jd₀ + t / 86400)

    # Obtain the position vector in the PEF reference frame to compute the gravity.
    r_pef = D_pef_tod * r_tod

    # Compute the gravity in the TOD reference frame.
    #
    # Since we are testing the J2 osculating propagator, we should use only the 2nd degree
    # and order 0 when computing the gravity. In this case, only the term related to the
    # constant J2 will be used.
    g_pef = compute_g(egm08, r_pef, 2, 0)
    g_tod = D_pef_tod' * g_pef

    # Compute the time derivative of the state vector.
    du = vcat(v_tod, g_tod)

    return du
end

############################################################################################
#                                           Run
############################################################################################

create_j2_osc_validation_tests()
