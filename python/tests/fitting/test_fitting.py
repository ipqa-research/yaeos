import yaeos
import yaeos.fitting
import yaeos.fitting.model_setters


def test_fit_co2_c6_pxy(data_co2_c6_pxy):
    data, tc, pc, w = data_co2_c6_pxy

    mr0 = yaeos.QMR(kij=[[0, 0], [0, 0]], lij=[[0, 0], [0, 0]])
    model = yaeos.PengRobinson76(tc, pc, w, mixrule=mr0)

    fit_kij = yaeos.fitting.model_setters.fit_kij_lij
    args = (model, True, False)

    problem = yaeos.fitting.BinaryFitter(
        model_setter=fit_kij, model_setter_args=args, data=data
    )

    x0 = [0.0]

    problem.fit(x0, bounds=None)

    assert abs(problem.solution.x[0] - 0.11) < 1e-1
