import glob

from pathlib import Path

from openalea.hydroroot.build_architecture import Config, process_plant


def test_process_plant():
    p = Path("data")

    cfg = Config(
        path = p.name,
        list_fn = ["plant1"]
    )

    directory = glob.glob(cfg.path)[0]
    prefix_fn = cfg.list_fn[0]
    d_r = process_plant(directory, prefix_fn, cfg)

    assert d_r['distance_from_base_(mm)'][0] == 10
    assert d_r['distance_from_base_(mm)'][1] == 20
    assert d_r['distance_from_base_(mm)'][2] == 80
    assert d_r['distance_from_base_(mm)'][3] == 90
    assert d_r['distance_from_base_(mm)'][4] == 158