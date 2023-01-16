within Workspace;

model ClosedSteamCycleTest
  extends .Modelon.Icons.Experiment;
  replaceable package Medium = .Modelon.Media.PreDefined.TwoPhase.WaterIF97 constrainedby
    .Modelon.Media.Interfaces.TwoPhaseMedium annotation(choicesAllMatching);
  parameter .Modelica.Units.SI.MassFlowRate steamHPflow=71 "Nominal HP steam flow rate";
  parameter .Modelica.Units.SI.AbsolutePressure p_steamHP=130e5 "Nominal HP steam pressure";
  parameter .Modelica.Units.SI.AbsolutePressure steamHPNomPressure=129.6e5 "Nominal HP steam pressure";

  parameter .Modelica.Units.SI.AbsolutePressure steamLPNomPressure=0.5e5 "Nominal IP steam pressure";

  parameter .Modelica.Units.SI.SpecificEnthalpy hl_HP=Medium.bubbleEnthalpy_pX(p_steamHP)
    "Start value of HP bubble point enthalpy";
  parameter .Modelica.Units.SI.SpecificEnthalpy hv_HP=Medium.dewEnthalpy_pX(p_steamHP)
    "Start value of HP dew point enthalpy";
  parameter .Modelica.Units.SI.Temperature T_HP=Medium.saturationTemperature_pX(p_steamHP)
    "Start value of HP saturation temperature";
  .Modelica.Units.SI.AbsolutePressure visualizer_p[5]={volume_afterDrum.p,condenser.p,condenser.p,pump_condenser.drain.p,
      volume_afterDrum.p};
  .Modelica.Units.SI.SpecificEnthalpy visualizer_h[5]={volume_afterDrum.h,condenser.summary.h_in[1],condenser.summary.h_out,
      pump_condenser.h,volume_afterDrum.h};
public
  .ThermalPower.TwoPhase.FlowResistances.HeightDiff downcomer(
    height=-10,
    useLoss=false,
    redeclare package Medium = Medium) annotation (Placement(transformation(
        origin={-5,-50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  .ThermalPower.TwoPhase.Drum drum(
    pstart=p_steamHP,
    hvstart=hv_HP,
    hlstart=hl_HP,
    ystart=0,
    Cm=0,
    DrumOrientation=1,
    rint=1,
    rext=1.1,
    L=2,
    Tmstart=600,
    initOpt=.Modelon.ThermoFluid.Choices.InitOptions.initialValues,
    redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-210,0},{-170,40}}, rotation=
            0)));
  .ThermalPower.TwoPhase.Valves.ValveSteam steamValve(
    CvData=.ThermalPower.SubComponents.Internal.Choices.FlowCoefficients.OpPoint,
    m_flow_nom=steamHPflow,
    p_nom=p_steamHP,
    T_nom=T_HP + 1,
    CheckValve=true,
    dp_nom=200000,
    redeclare package Medium = Medium)
                   annotation (Placement(transformation(extent={{-110,50},{-90,
            70}}, rotation=0)));

  .ThermalPower.Thermal.Sources.HeatFlowSource boiler(N=1) annotation (Placement(
        transformation(
        origin={-153,-42.5},
        extent={{9.5,-10},{-9.5,10}},
        rotation=90)));
  .ThermalPower.ControllersAndSensors.UniversalSensor levelSP1(outValue=0) annotation (
      Placement(transformation(extent={{-51,-51},{-61,-41}}, rotation=0)));
  .ThermalPower.ControllersAndSensors.UniversalSensor level1(outValue=drum.y) annotation (
      Placement(transformation(extent={{-62,-66},{-72,-56}}, rotation=0)));
  .ThermalPower.TwoPhase.TurboMachinery.Turbines.SteamTurbineStodola turbine(
    pstart=steamLPNomPressure,
    use_partialArc=false,
    p1_nom=steamHPNomPressure,
    m_flow_nom=steamHPflow,
    pstartin=steamHPNomPressure,
    p2_nom=steamLPNomPressure,
    hstartin=hl_HP,
    hstartout=2500e3,
    useNominalPoint=true,
    redeclare package Medium = Medium)
                          annotation (Placement(transformation(extent={{-61,36},
            {-31,66}},
          rotation=0)));
  .ThermalPower.TwoPhase.Volumes.MixVolume volume_preTurbine(
    N_drain=1,
    pstart=steamHPNomPressure,
    hstart=hv_HP,
    V_tot=1,
    N_feed=1,
    redeclare package Medium = Medium)
              annotation (Placement(transformation(extent={{-83,53},{-69,67}})));
  .Modelon.ThermoFluid.Electrical.SimpleGenerator generator
    annotation (Placement(transformation(extent={{-11,34},{23,68}})));
  .Modelica.Blocks.Sources.Ramp boilingHeat(
    duration=10,
    offset=100e6,
    height=20e6,
    startTime=500) annotation (Placement(transformation(extent={{-123,-36},{
            -139,-20}},
                   rotation=0)));
  .ThermalPower.TwoPhase.Condensers.Condenser condenser(
    redeclare package Medium = Medium,
    diameter=3,
    pstart=steamLPNomPressure,
    redeclare package CoolMedium = Medium,
    N_tubes=5000,
    ystart=1,
    V_hotwell=0,
    length=7,
    states=.ThermalPower.SubComponents.Internal.Choices.ThermoStates.UM)
              annotation (Placement(transformation(extent={{-20,-20},{10,11}},
          rotation=0)));
  .ThermalPower.TwoPhase.SourcesAndSinks.MassFlowBoundary coolFlow(
    m_flow0=7000,
    h0=Medium.specificEnthalpy(Medium.setState_pTX(
        1.01e5,
        273 + 15,
        Medium.fixedComposition)),
    redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{40,0},{20,20}}, rotation=0)));
  .ThermalPower.TwoPhase.SourcesAndSinks.PressureBoundary_h coolSink(N_ports=1, redeclare package
              Medium = Medium) annotation (Placement(transformation(extent={{40,-30},
            {20,-10}},      rotation=0)));
  .Modelica.Blocks.Continuous.LimPID drumLevelController(
    controllerType=.Modelica.Blocks.Types.SimpleController.PI,
    yMax=50000,
    yMin=0,
    initType=.Modelica.Blocks.Types.Init.InitialOutput,
    y_start=1500,
    k=500,
    Ti=25,
    Td=0.1) annotation (Placement(transformation(extent={{-67,-54},{-83,-38}})));
  .ThermalPower.TwoPhase.TurboMachinery.Pumps.PumpPosDispl pump_condenser(
    q_nom=steamHPflow/1000,
    m_flow_nom=steamHPflow,
    pin_start=steamLPNomPressure + 1e5,
    pout_start=steamHPNomPressure + 10e5,
    V=0.1,
    T_inertia=0.1,
    redeclare package Medium = Medium,
    redeclare package SatMedium = Medium,
    use_in_Np=false)
    annotation (Placement(transformation(extent={{-89,-91},{-109,-71}})));
  .ThermalPower.TwoPhase.Sensors.MultiData multiData(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-120,-85},{-140,-65}})));
  .ThermalPower.TwoPhase.Sensors.MultiData multiData1(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-38,-91},{-58,-71}})));
  .ThermalPower.TwoPhase.Sensors.MultiData multiData2(redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-140,70},{-120,50}})));
  .ThermalPower.TwoPhase.Sensors.MultiData multiData3(redeclare package Medium = Medium)
                                        annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={-16,20})));
  .ThermalPower.Visualizers.MultiDisplayVis_phTmdot multiDisplayVis_phTmdot
    annotation (Placement(transformation(extent={{-144.0,36.0},{-116.0,64.0}},rotation = 0.0,origin = {0.0,0.0})));
  .ThermalPower.Visualizers.MultiDisplayVis_phTmdot multiDisplayVis_phTmdot1
    annotation (Placement(transformation(extent={{-49,13},{-21,41}})));
  .ThermalPower.Visualizers.MultiDisplayVis_phTmdot multiDisplayVis_phTmdot2
    annotation (Placement(transformation(extent={{-144,-78},{-116,-50}})));
  .ThermalPower.Visualizers.MultiDisplayVis_phTmdot multiDisplayVis_phTmdot3
    annotation (Placement(transformation(extent={{-62,-84},{-34,-56}})));
  .ThermalPower.TwoPhase.Sensors.MultiData multiData4(redeclare package Medium = Medium)
                                        annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=270,
        origin={-176,-3})));
  .ThermalPower.Visualizers.MultiDisplayVis_phTmdot multiDisplayVis_phTmdot4
    annotation (Placement(transformation(extent={{-169,-9},{-141,19}})));
  .ThermalPower.TwoPhase.TurboMachinery.Pumps.PumpPosDispl pump_riser(
    m_flow_nom=steamHPflow,
    pin_start=steamLPNomPressure + 1e5,
    pout_start=steamHPNomPressure + 10e5,
    q_nom=10*steamHPflow/1000,
    V=0.1,
    redeclare package Medium = Medium,
    redeclare package SatMedium = Medium,
    use_in_n=false,
    use_in_Np=false)
    annotation (Placement(transformation(extent={{-198,-70},{-178,-50}})));

  .ThermalPower.TwoPhase.FlowResistances.OrificeLiquid riser_orifice(
    CvData=.ThermalPower.SubComponents.Internal.Choices.FlowCoefficients.OpPoint,
    m_flow_nom=steamHPflow,
    d_nom=200,
    dp_nom=1000000,
    redeclare package Medium = Medium)
                    annotation (Placement(transformation(
        extent={{-6.5,-6.5},{6.5,6.5}},
        rotation=90,
        origin={-175.5,-25.5})));

  .ThermalPower.TwoPhase.Volumes.Header riser(
    N_drain=1,
    pstart=steamHPNomPressure + 10e5,
    hstart=(hv_HP + hl_HP)/2,
    V_tot=1,
    redeclare package Medium = Medium,
    N_feed=1)
             annotation (Placement(transformation(
        extent={{-7.5,8},{7.5,-8}},
        rotation=90,
        origin={-176,-42.5})));
  .ThermalPower.TwoPhase.Sensors.Pressure pSteam(redeclare package Medium = Medium, useSIunit=
       true) annotation (Placement(transformation(extent={{-147,73},{-134,85}},
          rotation=0)));
  .ThermalPower.ControllersAndSensors.UniversalSensor pSP(outValue=p_steamHP) annotation (
      Placement(transformation(extent={{-134,84},{-124,94}}, rotation=0)));
  .ThermalPower.ControllersAndSensors.LimPI pPI(
    k=-1e-6,
    steadyStateInit=false,
    yInit=1,
    Ti=2,
    yMin=0.01) annotation (Placement(transformation(extent={{-118,81},{-102,97}},
          rotation=0)));
  .ThermalPower.Visualizers.FourValueLegend fourValueLegend
    annotation (Placement(transformation(extent={{12,-85},{49,-52}})));

  // Masses
protected
  .Modelica.Units.SI.Mass m_drum=drum.Mv + drum.Ml;
  .Modelica.Units.SI.Mass m_volume_preTurbine=volume_preTurbine.M;
  .Modelica.Units.SI.Mass m_condenser=condenser.summary.m;
  .Modelica.Units.SI.Mass m_riser=riser.M;

  .Modelica.Units.SI.Mass m_tot=m_drum + m_volume_preTurbine + m_condenser + m_riser;
public
  .ThermalPower.TwoPhase.Volumes.Header volume_afterDrum(
    N_drain=2,
    pstart=steamHPNomPressure,
    hstart=hv_HP,
    V_tot=1,
    N_feed=1,
    redeclare package Medium = Medium)
    annotation (Placement(transformation(extent={{-172,50},{-152,70}})));
  .Modelica.Blocks.Sources.Ramp superHeat(
    duration=10,
    height=0,
    startTime=0,
    offset=20e6) annotation (Placement(transformation(extent={{-204,79},{-188,
            95}}, rotation=0)));
  .ThermalPower.Thermal.Sources.HeatFlowSource superHeater(N=1) annotation (Placement(
        transformation(
        origin={-162,85},
        extent={{10,-10},{-10,10}},
        rotation=180)));
  .ThermalPower.TwoPhase.FlowResistances.OrificeLiquid riser_orifice1(
    CvData=.ThermalPower.SubComponents.Internal.Choices.FlowCoefficients.OpPoint,
    m_flow_nom=steamHPflow,
    dp_nom=100000,
    d_nom=10,
    redeclare package Medium = Medium)
              annotation (Placement(transformation(
        extent={{-5.5,-6},{5.5,6}},
        rotation=90,
        origin={-179,50.5})));

  inner .ThermalPower.System_TPL system_TPL
    annotation (Placement(transformation(extent={{-230,40},{-210,60}})));
  .ThermalPower.Visualizers.PH_water pH_water(x=visualizer_h, y=visualizer_p)
    annotation (Placement(transformation(extent={{45,-41},{201,115}})));
equation
  connect(volume_preTurbine.drain[1], turbine.feed) annotation (Line(
      points={{-70.4,60},{-58,60}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(generator.flange_a, turbine.shaft_b) annotation (Line(
      points={{-7.6,51},{-32.5,51}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(boilingHeat.y, boiler.power) annotation (Line(
      points={{-139.8,-28},{-147,-28},{-147,-33.95}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(coolFlow.port, condenser.feed_cool) annotation (Line(
      points={{21,10},{11.8,10},{11.8,-0.16}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(condenser.drain_cool, coolSink.port[1]) annotation (Line(
      points={{11.8,-6.36},{11.8,-20},{21,-20}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(downcomer.feed, condenser.drain) annotation (Line(
      points={{-5,-40},{-5,-13.8}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(level1.y, drumLevelController.u_m) annotation (Line(
      points={{-72,-61},{-75,-61},{-75,-55.6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(levelSP1.y, drumLevelController.u_s) annotation (Line(
      points={{-61,-46},{-65.4,-46}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(drumLevelController.y, pump_condenser.in_n) annotation (Line(
      points={{-83.8,-46},{-90,-46},{-90,-73},{-91.8,-73}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pump_condenser.drain, multiData.port_a) annotation (Line(
      points={{-107,-75},{-124,-75}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(multiData1.port_b, pump_condenser.feed) annotation (Line(
      points={{-54,-81},{-91,-81}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(downcomer.drain, multiData1.port_a) annotation (Line(
      points={{-5,-60},{-5,-81},{-42,-81}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(steamValve.feed, multiData2.port_b) annotation (Line(
      points={{-110,60},{-124,60}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(multiData3.port_a, turbine.drain) annotation (Line(
      points={{-16,26},{-16,40.5},{-34,40.5}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(multiData2.u, multiDisplayVis_phTmdot.y) annotation (Line(
      points={{-130,60},{-130,50}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(multiData3.u, multiDisplayVis_phTmdot1.y) annotation (Line(
      points={{-16,20},{-35,20},{-35,27}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(multiDisplayVis_phTmdot2.y, multiData.u) annotation (Line(
      points={{-130,-64},{-130,-75}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(multiDisplayVis_phTmdot3.y, multiData1.u) annotation (Line(
      points={{-48,-70},{-48,-81}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(drum.riser, multiData4.port_b) annotation (Line(
      points={{-176,10},{-176,3}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(multiDisplayVis_phTmdot4.y, multiData4.u) annotation (Line(
      points={{-155,5},{-155,-3},{-176,-3}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(multiData.port_b, drum.feed) annotation (Line(
      points={{-136,-75},{-208,-75},{-208,19.6}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(pump_riser.feed, drum.downcomer) annotation (Line(
      points={{-196,-60},{-203,-60},{-203,10},{-204,10}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(multiData4.port_a, riser_orifice.drain) annotation (Line(
      points={{-176,-9},{-176,-19},{-175.5,-19}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(riser.drain[1], riser_orifice.feed) annotation (Line(
      points={{-176,-36.5},{-176,-36},{-175.5,-36},{-175.5,-32}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(pSteam.p, pPI.u_m) annotation (Line(
      points={{-134.65,79},{-110,79},{-110,79.4}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pPI.y, steamValve.command) annotation (Line(
      points={{-101.2,89},{-100,89},{-100,65.6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(volume_preTurbine.feed[1], steamValve.drain) annotation (Line(
      points={{-81.6,60},{-90,60}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(volume_afterDrum.wall, superHeater.port[1]) annotation (Line(
      points={{-162,69},{-162,75}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(pSP.y, pPI.u_s) annotation (Line(
      points={{-124,89},{-119.6,89}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(drum.drain, riser_orifice1.feed) annotation (Line(
      points={{-178,31.6},{-178,45},{-179,45}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(riser_orifice1.drain, volume_afterDrum.feed[1]) annotation (Line(
      points={{-179,56},{-179,60},{-170,60}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(superHeat.y, superHeater.power) annotation (Line(
      points={{-187.2,87},{-182,87},{-182,91},{-171,91}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(volume_afterDrum.drain[1], pSteam.port) annotation (Line(
      points={{-154,59.75},{-154,73},{-140.5,73}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(multiData3.port_b, condenser.feed) annotation (Line(
      points={{-16,14},{-16,9},{-15.95,9},{-15.95,3.87}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(riser.feed[1], pump_riser.drain) annotation (Line(points={{-176,-48.5},
          {-176,-54},{-180,-54}}, color={0,0,255}));
  connect(riser.wall, boiler.port[1]) annotation (Line(points={{-168.8,-42.5},{
          -163,-42.5},{-163,-42.5}},   color={191,0,0}));
  connect(volume_afterDrum.drain[2], multiData2.port_a) annotation (Line(points={{-154,
          60.25},{-154,60},{-136,60}},      color={0,0,255}));
  annotation (
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-250,-100},{200,100}},
        grid={1,1})),
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        grid={1,1})),
    Documentation(info="<html>
<h4>Description</h4>
<p>This example demonstrates a closed steam cycle. The drum level and pressure is controlled by adjusting the feed water pump respective the high pressure steam valve.</p>
<p>There are some transients in the beginning of the simulation due to the system is not started in steady state. </p>
<p>At t=500 there is an increase of the heat load to the riser tubes. This will increase the pressure initially before the pressure controller controls the pressure back to the set-point by increasing the high pressure steam valve opening which result in a higher steam flow and an increase of the produced power.</p>
<h4>Experiment setup</h4>
<p>Simulate for 1000s with a tolerance of 1e-5.</p>
<h4>Output</h4>
<p>Variables of interest are:</p>
<ul>
<li>produced power, variable power in component generator</li>
<li>steam flow rate in the turbine, variable summary.m_flow in component turbine</li>
<li>pressure, temperature and level in the drum, variables p,T and y in component drum</li>
<li>Condenser pressure and heat, variables p and Q_cool in component condenser</li>
<li>heat to boiler, variable power in component boiler</li>
</ul>
</html>", revisions="<html>
Copyright &copy; 2004-2022, MODELON AB <br /> The use of this software component is regulated by the licensing conditions for Modelon Libraries. <br />This copyright notice must, unaltered, accompany all components that are derived from, copied from, <br />or by other means have their origin from any Modelon Library.
</html>"),
    experiment(StopTime=1000, Tolerance=1e-006));
end ClosedSteamCycleTest;
