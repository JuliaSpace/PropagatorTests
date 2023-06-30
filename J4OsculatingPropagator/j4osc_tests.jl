# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#                                    Author : Ronan Arraes Jardim Chagas <ronisbr@gmail.com>
#                                      Date : 2023-06-30T13:30:00 UTC-3
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Create a test set for the J4 osculating orbit propagator using a numerical algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

using Dates
using DifferentialEquations
using PrettyTables
using SatelliteToolboxGravityModels
using SatelliteToolboxTransformations
using StaticArrays

"""
    create_j4osc_validation_tests()

Create the validation tests for the J4 osculating orbit propagator.
"""
function create_j4osc_validation_tests()
    # Load gravity model
    # ======================================================================================

    # We are using the constants in the EGM-2008 model
    egm2008 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM2008))

    # Inputs
    # ======================================================================================

    # This state vector is obtained using the J4 osculating propagator with `t = 0s` and the
    # following mean elements:
    #
    #   - Epoch           : 2023-01-01T00:00:00
    #   - Semi-major axis : 7190.982    km
    #   - Eccentricity    :    0.001111
    #   - Inclination     :   98.405    °
    #   - RAAN            :   90.0      °
    #   - Arg. of Perigee :  200.0      °
    #   - True Anomaly    :   45.0      °

    jd₀ = DateTime("2023-01-01") |> datetime2julian
    r₀  = @SVector [-952882.6807035431, -3.03843825141778e6, -6.444903144699051e6]
    v₀  = @SVector [-460.11262511481317, 6745.426195633091, -3115.9662885215403]

    # Build the ODE problem
    # ======================================================================================

    p = (; jd₀, egm2008)
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
    jd₀     = p.jd₀
    egm2008 = p.egm2008

    # Unpack the state.
    r_tod = SVector(u[1], u[2], u[3])
    v_tod = SVector(u[4], u[5], u[6])

    # Compute the matrix to convert the TOD reference frame to the PEF reference frame.
    D_pef_tod = r_eci_to_ecef(TOD(), PEF(), jd₀ + t / 86400)

    # Obtain the position vector in the PEF reference frame to compute the gravity.
    r_pef = D_pef_tod * r_tod

    # Compute the gravity in the TOD reference frame.
    #
    # Since we are testing the J4 osculating propagator, we should use only the 2nd degree
    # and order 0 when computing the gravity. In this case, only the term related to the
    # constant J4 will be used.
    g_pef = GravityModels.gravitational_acceleration(
        egm2008,
        r_pef;
        max_degree = 4,
        max_order  = 0
    )
    g_tod = D_pef_tod' * g_pef

    # Compute the time derivative of the state vector.
    du = vcat(v_tod, g_tod)

    return du
end

############################################################################################
#                                           Run
############################################################################################

create_j4osc_validation_tests()
