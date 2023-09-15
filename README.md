# End-to-end-process-monitoring
Simulation used to demonstrate end to end process monitoring, in preparation for IFAC World Congress 2023

To do:
* Add griddedInterpolant function
* Implement control charts as alternate monitoring method
* Implement two different advisories
* Separate "r" into a supervisory control signal "r" and a maintenance signal "q"
* Split "r.components" into "r.Sensors" and "r.Valves". This will have the type built-in directly, and can be used to generalize the regulatory control module
* Generalize the supervisory control module:
  * Define "r.Shut.FW.Position", then use "r.Valve.FW.Position = r.(r.regime).FW.Position", which will move the defintion of the regimes to the parameter deifnition phase
  * Let "r.setpoints.C(end+1) = r.(r.regime).setpoints.C" (same as above)
* Generalize the monitoring module and improve the overall monitoring methodology
