import os
import knoblist

knobNames = knoblist.GENIEMultisigmaKnobNames

for name in knobNames:

  IsEnabled = 'true'
  if 'MFP' in name or 'Fr' in name:
    IsEnabled = 'false'

  output = '''- parameterName: "%s"
  isEnabled: %s
  dialSetDefinitions:
    - dialType: Spline
      minimumSplineResponse: 0
      dialLeafName: "%s"
      applyCondition: "[reco_is_data]==0"
'''%(name, IsEnabled, name)
  print(output)

knobNames = knoblist.GENIEMorphKnobNames

for name in knobNames:
  output = '''- parameterName: "%s"
  isEnabled: true
  dialSetDefinitions:
    - dialType: Spline
      minimumSplineResponse: 0
      dialLeafName: "%s"
      useMirrorDial: true
      mirrorLowEdge: -1
      mirrorHighEdge: 1
      applyCondition: "[reco_is_data]==0"
'''%(name, name)
  print(output)
