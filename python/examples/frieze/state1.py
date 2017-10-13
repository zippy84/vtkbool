# state file generated using paraview version 5.4.1

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [894, 612]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [-49.29166666666667, 68.25, 0.541665]
renderView1.StereoType = 0
renderView1.CameraPosition = [-69.3183238895627, -145.49349251302502, 261.94215377564524]
renderView1.CameraFocalPoint = [5.316607220996967, 10.100912301201648, -17.265806958328785]
renderView1.CameraViewUp = [0.5563317767903634, 0.6537008872804081, 0.5130010761217609]
renderView1.CameraParallelScale = 10.436136517218301
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.32, 0.34, 0.43]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
legacyVTKReader1 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test5.vtk'])

# create a new 'Legacy VTK Reader'
legacyVTKReader3 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test4.vtk'])

# create a new 'Transform'
transform1 = Transform(Input=legacyVTKReader1)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [-98.58333333333333, 0.0, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, -90.0]

# create a new 'Legacy VTK Reader'
legacyVTKReader4 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test7.vtk'])

# create a new 'Transform'
transform4 = Transform(Input=legacyVTKReader4)
transform4.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform4.Transform.Translate = [-98.5833333, 136.5, 0.0]
transform4.Transform.Rotate = [0.0, 0.0, 180.0]

# create a new 'Legacy VTK Reader'
legacyVTKReader2 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test8.vtk'])

# create a new 'Transform'
transform2 = Transform(Input=legacyVTKReader2)
transform2.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform2.Transform.Translate = [-95.3333333, 179.8333, 0.0]
transform2.Transform.Rotate = [0.0, 0.0, -90.0]

# create a new 'Legacy VTK Reader'
test14vtk = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test14.vtk'])

# create a new 'Transform'
transform5 = Transform(Input=test14vtk)
transform5.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform5.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Legacy VTK Reader'
legacyVTKReader5 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test6.vtk'])

# create a new 'Transform'
transform3 = Transform(Input=legacyVTKReader5)
transform3.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform3.Transform.Translate = [-98.58333333333333, 136.5, 0.0]
transform3.Transform.Rotate = [0.0, 0.0, -90.0]

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get color transfer function/color map for 'RegionId'
regionIdLUT = GetColorTransferFunction('RegionId')
regionIdLUT.RGBPoints = [0.0, 0.231373, 0.298039, 0.752941, 1.0, 0.865003, 0.865003, 0.865003, 2.0, 0.705882, 0.0156863, 0.14902]
regionIdLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'RegionId'
regionIdPWF = GetOpacityTransferFunction('RegionId')
regionIdPWF.Points = [0.0, 0.0, 0.5, 0.0, 2.0, 1.0, 0.5, 0.0]
regionIdPWF.ScalarRangeInitialized = 1

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from legacyVTKReader3
legacyVTKReader3Display = Show(legacyVTKReader3, renderView1)
# trace defaults for the display properties.
legacyVTKReader3Display.Representation = 'Surface With Edges'
legacyVTKReader3Display.ColorArrayName = ['POINTS', 'RegionId']
legacyVTKReader3Display.LookupTable = regionIdLUT
legacyVTKReader3Display.EdgeColor = [0.0, 0.0, 0.0]
legacyVTKReader3Display.OSPRayScaleArray = 'RegionId'
legacyVTKReader3Display.OSPRayScaleFunction = 'PiecewiseFunction'
legacyVTKReader3Display.SelectOrientationVectors = 'None'
legacyVTKReader3Display.ScaleFactor = 9.966666700000001
legacyVTKReader3Display.SelectScaleArray = 'RegionId'
legacyVTKReader3Display.GlyphType = 'Arrow'
legacyVTKReader3Display.GlyphTableIndexArray = 'RegionId'
legacyVTKReader3Display.DataAxesGrid = 'GridAxesRepresentation'
legacyVTKReader3Display.PolarAxes = 'PolarAxesRepresentation'
legacyVTKReader3Display.GaussianRadius = 4.9833333500000006
legacyVTKReader3Display.SetScaleArray = ['POINTS', 'RegionId']
legacyVTKReader3Display.ScaleTransferFunction = 'PiecewiseFunction'
legacyVTKReader3Display.OpacityArray = ['POINTS', 'RegionId']
legacyVTKReader3Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
legacyVTKReader3Display.SetScalarBarVisibility(renderView1, True)

# show data from transform1
transform1Display = Show(transform1, renderView1)
# trace defaults for the display properties.
transform1Display.Representation = 'Surface With Edges'
transform1Display.ColorArrayName = ['POINTS', 'RegionId']
transform1Display.LookupTable = regionIdLUT
transform1Display.EdgeColor = [0.0, 0.0, 0.0]
transform1Display.OSPRayScaleArray = 'RegionId'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 9.3624967
transform1Display.SelectScaleArray = 'RegionId'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'RegionId'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.GaussianRadius = 4.68124835
transform1Display.SetScaleArray = ['POINTS', 'RegionId']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', 'RegionId']
transform1Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform1Display.SetScalarBarVisibility(renderView1, True)

# show data from transform3
transform3Display = Show(transform3, renderView1)
# trace defaults for the display properties.
transform3Display.Representation = 'Surface With Edges'
transform3Display.ColorArrayName = ['POINTS', 'RegionId']
transform3Display.LookupTable = regionIdLUT
transform3Display.EdgeColor = [0.0, 0.0, 0.0]
transform3Display.OSPRayScaleArray = 'RegionId'
transform3Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform3Display.SelectOrientationVectors = 'None'
transform3Display.ScaleFactor = 4.3874967
transform3Display.SelectScaleArray = 'RegionId'
transform3Display.GlyphType = 'Arrow'
transform3Display.GlyphTableIndexArray = 'RegionId'
transform3Display.DataAxesGrid = 'GridAxesRepresentation'
transform3Display.PolarAxes = 'PolarAxesRepresentation'
transform3Display.GaussianRadius = 2.19374835
transform3Display.SetScaleArray = ['POINTS', 'RegionId']
transform3Display.ScaleTransferFunction = 'PiecewiseFunction'
transform3Display.OpacityArray = ['POINTS', 'RegionId']
transform3Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform3Display.SetScalarBarVisibility(renderView1, True)

# show data from transform4
transform4Display = Show(transform4, renderView1)
# trace defaults for the display properties.
transform4Display.Representation = 'Surface With Edges'
transform4Display.ColorArrayName = ['POINTS', 'RegionId']
transform4Display.LookupTable = regionIdLUT
transform4Display.EdgeColor = [0.0, 0.0, 0.0]
transform4Display.OSPRayScaleArray = 'RegionId'
transform4Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform4Display.SelectOrientationVectors = 'None'
transform4Display.ScaleFactor = 0.47916670000000006
transform4Display.SelectScaleArray = 'RegionId'
transform4Display.GlyphType = 'Arrow'
transform4Display.GlyphTableIndexArray = 'RegionId'
transform4Display.DataAxesGrid = 'GridAxesRepresentation'
transform4Display.PolarAxes = 'PolarAxesRepresentation'
transform4Display.GaussianRadius = 0.23958335000000003
transform4Display.SetScaleArray = ['POINTS', 'RegionId']
transform4Display.ScaleTransferFunction = 'PiecewiseFunction'
transform4Display.OpacityArray = ['POINTS', 'RegionId']
transform4Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform4Display.SetScalarBarVisibility(renderView1, True)

# show data from transform2
transform2Display = Show(transform2, renderView1)
# trace defaults for the display properties.
transform2Display.Representation = 'Surface With Edges'
transform2Display.ColorArrayName = ['POINTS', 'RegionId']
transform2Display.LookupTable = regionIdLUT
transform2Display.EdgeColor = [0.0, 0.0, 0.0]
transform2Display.OSPRayScaleArray = 'RegionId'
transform2Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform2Display.SelectOrientationVectors = 'None'
transform2Display.ScaleFactor = 3.466663
transform2Display.SelectScaleArray = 'RegionId'
transform2Display.GlyphType = 'Arrow'
transform2Display.GlyphTableIndexArray = 'RegionId'
transform2Display.DataAxesGrid = 'GridAxesRepresentation'
transform2Display.PolarAxes = 'PolarAxesRepresentation'
transform2Display.GaussianRadius = 1.7333315
transform2Display.SetScaleArray = ['POINTS', 'RegionId']
transform2Display.ScaleTransferFunction = 'PiecewiseFunction'
transform2Display.OpacityArray = ['POINTS', 'RegionId']
transform2Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform2Display.SetScalarBarVisibility(renderView1, True)

# show data from transform5
transform5Display = Show(transform5, renderView1)
# trace defaults for the display properties.
transform5Display.Representation = 'Surface With Edges'
transform5Display.ColorArrayName = ['POINTS', 'RegionId']
transform5Display.LookupTable = regionIdLUT
transform5Display.EdgeColor = [0.0, 0.0, 0.0]
transform5Display.OSPRayScaleArray = 'RegionId'
transform5Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform5Display.SelectOrientationVectors = 'None'
transform5Display.ScaleFactor = 1.2999967000000001
transform5Display.SelectScaleArray = 'RegionId'
transform5Display.GlyphType = 'Arrow'
transform5Display.GlyphTableIndexArray = 'RegionId'
transform5Display.DataAxesGrid = 'GridAxesRepresentation'
transform5Display.PolarAxes = 'PolarAxesRepresentation'
transform5Display.GaussianRadius = 0.6499983500000001
transform5Display.SetScaleArray = ['POINTS', 'RegionId']
transform5Display.ScaleTransferFunction = 'PiecewiseFunction'
transform5Display.OpacityArray = ['POINTS', 'RegionId']
transform5Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform5Display.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for regionIdLUT in view renderView1
regionIdLUTColorBar = GetScalarBar(regionIdLUT, renderView1)
regionIdLUTColorBar.Title = 'RegionId'
regionIdLUTColorBar.ComponentTitle = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(transform5)
# ----------------------------------------------------------------
