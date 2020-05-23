# chemhelp

Started out as common functionality I keep re-programming, now it is a collection of interfaces between RDKit and quantum chemistry

For example, want to optimize a molecule with PM3 in GAMESS?

    from chemhelp import gamess
    from chemhelp import cheminfo
    
    methane = """ SDF CONTENT """

    header = """ $basis gbasis=pm3 $end
     $contrl runtyp=optimize icharg={:} $end
     $statpt opttol=0.0005 nstep=300 projct=.F. $end """

    molobj = cheminfo.sdfstr_to_molobj(methane)
    stdout, stderr = gamess.calculate(molobj, header)
    properties = gamess.read_properties_coordinates(stdout)
    atoms = properties["atoms"]
    coord = properties["coord"]
    energy = properties["h"]

with consistent interfaces to MOPAC, GAMESS, Gaussian and MNDO

