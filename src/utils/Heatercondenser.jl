export OHPConfiguration

function OHPConfiguration(configure_type::String,power::Real,Tc::Real,Δx::Real;hc::Real=2300.0,hc2ratio=1/30)

    if configure_type == "ASETS-II OHP 1 LARGE HEATER"
        total_heater_area = 2.0inches*2.0inches;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.5inches,1.0inches,1.5*Δx)
        Tfe = RigidTransform((0.7inches,-0.0),0.0)
        Tfe(eb1)

        eb2 = Rectangle(0.5inches,1.0inches,1.5*Δx)
        Tfe = RigidTransform((-0.7inches,-0.0),0.0)
        Tfe(eb2)

        cb1 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx) # 0.02916 = 0.0648*0.9/2 
        Tfc = RigidTransform((-2.45inches,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx)
        Tfc = RigidTransform((2.45inches,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        eparams2 = PrescribedHeatFluxRegion(qe,eb2);
        cparams1 = PrescribedHeatModelRegion(hc,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1,eparams2], [cparams1,cparams2]
    end

    # In one sided condenser case, the "adiabatic" side is not completely adiabatic in the code. 
    # Instead it is a fraction "hc2ratio" of the regular condenser htc. 
    # And it is a tunable parameter for now as it is a representaion of the insulation material:)

    if (configure_type == "ASETS-II OHP 2 LARGE HEATER") || (configure_type == "ASETS-II OHP 3 LARGE HEATER")
        total_heater_area = 2.0inches*2.0inches;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.5inches,1.0inches,1.5*Δx)
        Tfe = RigidTransform((0.7inches,-0.0),0.0)
        Tfe(eb1)

        eb2 = Rectangle(0.5inches,1.0inches,1.5*Δx)
        Tfe = RigidTransform((-0.7inches,-0.0),0.0)
        Tfe(eb2)

        cb1 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx) # 0.02916 = 0.0648*0.9/2 
        Tfc = RigidTransform((-2.45inches,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx)
        Tfc = RigidTransform((2.45inches,-0.0),0.0)
        Tfc(cb2)

        # cb3 = Rectangle(1.9inches,1.2inches ,1.5*Δx)
        # Tfc = RigidTransform((0.0,-0.0),0.0)
        # Tfc(cb3)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        eparams2 = PrescribedHeatFluxRegion(qe,eb2);
        cparams1 = PrescribedHeatModelRegion(hc*hc2ratio,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);
        # cparams3 = PrescribedHeatModelRegion(10.0,Tc,cb3);

    return [eparams1,eparams2], [cparams1,cparams2]
    end

    if configure_type == "ASETS-II OHP 1 SMALL HEATER"
        total_heater_area = 0.5inches*0.5inches;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.25inches,0.25inches,1.5*Δx)
        Tfe = RigidTransform((0.0inches,-0.0),0.0)
        Tfe(eb1)

        cb1 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx) # 0.02916 = 0.0648*0.9/2 
        Tfc = RigidTransform((-2.45inches,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx)
        Tfc = RigidTransform((2.45inches,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        cparams1 = PrescribedHeatModelRegion(hc,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1], [cparams1,cparams2]
    end

    if (configure_type == "ASETS-II OHP 2 SMALL HEATER") || (configure_type == "ASETS-II OHP 3 SMALL HEATER")
        total_heater_area = 0.5inches*0.5inches;
        qe = power/total_heater_area;

        eb1 = Rectangle(0.25inches,0.25inches,1.5*Δx)
        Tfe = RigidTransform((0.0inches,-0.0),0.0)
        Tfe(eb1)

        cb1 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx) # 0.02916 = 0.0648*0.9/2 
        Tfc = RigidTransform((-2.45inches,-0.0),0.0)
        Tfc(cb1)

        cb2 = Rectangle(0.5inches,0.0648*0.95/2 ,1.5*Δx)
        Tfc = RigidTransform((2.45inches,-0.0),0.0)
        Tfc(cb2)

        eparams1 = PrescribedHeatFluxRegion(qe,eb1);
        cparams1 = PrescribedHeatModelRegion(hc*hc2ratio,Tc,cb1);
        cparams2 = PrescribedHeatModelRegion(hc,Tc,cb2);

    return [eparams1], [cparams1,cparams2]
    end


    return "configuration not recognized"
end

