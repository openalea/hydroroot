import glob

from pathlib import Path

from openalea.hydroroot.build_architecture import Config, Results, process_plant


def test_process_plant():
    p = Path("data")
    res = Results()
    cfg = Config(
        path = p.name,
        list_fn = ["plant1"]
    )

    process_plant(cfg, res)

    assert res.d_archi[cfg.list_fn[0]]['distance_from_base_(mm)'][0] == 10
    assert res.d_archi[cfg.list_fn[0]]['distance_from_base_(mm)'][1] == 20
    assert res.d_archi[cfg.list_fn[0]]['distance_from_base_(mm)'][2] == 80
    assert res.d_archi[cfg.list_fn[0]]['distance_from_base_(mm)'][3] == 90
    assert res.d_archi[cfg.list_fn[0]]['distance_from_base_(mm)'][4] == 158