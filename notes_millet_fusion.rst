**Duplicates removed**

- millet.conductance.compute_k == hydroroot.conductance.cimpute_k
- openalea.hydroroot.millet.conductance.compute_K -> openalea.hydroroot.conductance.compute_K_from_Poiseuille not used
- openalea.hydroroot.millet.conductance.poiseuille == openalea.hydroroot.conductance.poiseuille used in above
- openalea.hydroroot.millet.conductance.fit_property_from_csv == openalea.hydroroot.conductance.fit_property_from_csv
- openalea.hydroroot.millet.conductance.fit_property == openalea.hydroroot.conductance.fit_property
- openalea.hydroroot.millet.conductance.fit_K == openalea.hydroroot.conductance.fit_K
- openalea.hydroroot.millet.length.fit_length == openalea.hydroroot.length.fit_length
- openalea.hydroroot.millet.length.fit_law == openalea.hydroroot.length.fit_law
- openalea.hydroroot.millet.radius.cont_radius == openalea.hydroroot.radius.cont_radius and none are used
- openalea.hydroroot.millet.radius.discont_radius == openalea.hydroroot.radius.discont_radius
- openalea.hydroroot.millet.radius.ordered_radius == openalea.hydroroot.radius.ordered_radius
- openalea.hydroroot.millet.radius.compute_length == openalea.hydroroot.radius.compute_length
- openalea.hydroroot.millet.radius.compute_relative_position == openalea.hydroroot.radius.compute_relative_position
- openalea.hydroroot.millet.radius.compute_surface == openalea.hydroroot.radius.compute_surface
- openalea.hydroroot.millet.radius.compute_volume == openalea.hydroroot.radius.compute_volume

**Moved**

- openalea.hydroroot.millet.radius.radius not used moved to hydroroot.radius


**Kept but duplicates with hydroroot**

- openalea.hydroroot.millet.conductance.fit_property_from_spline used in millet.law.add_soil that is used in example/millet/Water flux-Copy1.ipynb, very similar with hydroroot the only difference `dict(list(zip(keys, y_values)))` in hydroroot
- openalea.hydroroot.law.expovariate_law is a kind of duplicate of openalea.hydroroot.millet.law.expovariate_law + openalea.hydroroot.millet.law.read_data. None of them are used but read_data is used by openalea.hydroroot.millet.law.age_law that is not used
    - big difference millet law gave things in nb of vertices, hydroroot in meter, note that none of them are correct because of discretize function see below
- openalea.hydroroot.millet.law.discretize is not sctrickly a duplicate but return the same than openalea.hydroroot.law.discretize, and it is used only in openalea.hydroroot.millet.law.expovariate_law that is not used.
    - both descritize function seem not correct to be verified
- openalea.hydroroot.millet.law.compute_age_with_constant_growth_speed is a duplicate of openalea.hydroroot.millet.age.compute_age_with_constant_growth_speed and there are not used
- openalea.hydroroot.millet.law.radius_from_computed_diameters == openalea.hydroroot.radius.radius we should change the later by the former

***Remark:***

- openalea.hydroroot.conductance.fit_property_from_csv is deprecated and use pylab I think because used in wralea
- openalea.hydroroot.conductance.fit_property used in conductance.fit_K that is not used and use pylab
- plot function in src/openalea/hydroroot/millet/display_millet.py not duplicate but used millet.law and millet.turtle to account for time 
