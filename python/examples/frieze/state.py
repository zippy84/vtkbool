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
renderView1.CenterOfRotation = [2.4165000915527344, 1.875, 0.5416499972343445]
renderView1.StereoType = 0
renderView1.CameraPosition = [-29.979463243401923, -171.61649028674506, 73.62708490569474]
renderView1.CameraFocalPoint = [-92.7198238963452, -38.72104348302699, 9.488855898247008]
renderView1.CameraViewUp = [-0.07405037209849341, 0.40488478206314543, 0.9113642826256436]
renderView1.CameraParallelScale = 6.168850567681643
renderView1.CameraParallelProjection = 1
renderView1.Background = [0.32, 0.34, 0.43]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
legacyVTKReader1 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test10.vtk'])

# create a new 'Legacy VTK Reader'
legacyVTKReader3 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test1.vtk'])

# create a new 'Legacy VTK Reader'
legacyVTKReader4 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test12.vtk'])

# create a new 'Legacy VTK Reader'
legacyVTKReader7 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test9.vtk'])

# create a new 'Transform'
transform2 = Transform(Input=legacyVTKReader4)
transform2.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform2.Transform.Translate = [-102.9167, -21.6666, 0.0]
transform2.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Legacy VTK Reader'
legacyVTKReader2 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test3.vtk'])

# create a new 'Legacy VTK Reader'
legacyVTKReader5 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test2.vtk'])

# create a new 'Legacy VTK Reader'
legacyVTKReader6 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test11.vtk'])

# create a new 'Transform'
transform3 = Transform(Input=legacyVTKReader7)
transform3.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform3.Transform.Translate = [-39.0, -7.583333333333333, 0.0]
transform3.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Legacy VTK Reader'
legacyVTKReader8 = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test0.vtk'])

# create a new 'Transform'
transform7 = Transform(Input=legacyVTKReader1)
transform7.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform7.Transform.Translate = [-52.0, -20.5832, 0.0]

# create a new 'Transform'
transform4 = Transform(Input=legacyVTKReader6)
transform4.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform4.Transform.Translate = [-39.0, -7.5832, 0.0]
transform4.Transform.Rotate = [0.0, 0.0, 45.0]

# create a new 'Transform'
transform1 = Transform(Input=legacyVTKReader2)
transform1.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform1.Transform.Translate = [43.3333, 3.25, 0.0]
transform1.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Legacy VTK Reader'
test13vtk = LegacyVTKReader(FileNames=['/home/zippy/vtkbool_light/python/examples/frieze/test13.vtk'])

# create a new 'Transform'
transform8 = Transform(Input=test13vtk)
transform8.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform8.Transform.Translate = [-102.9167, -21.6666, 0.0]

# create a new 'Transform'
transform5 = Transform(Input=legacyVTKReader3)
transform5.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform5.Transform.Rotate = [0.0, 0.0, 90.0]

# create a new 'Transform'
transform6 = Transform(Input=legacyVTKReader5)
transform6.Transform = 'Transform'

# init the 'Transform' selected for 'Transform'
transform6.Transform.Translate = [43.333, 3.25, 0.0]

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

# show data from legacyVTKReader8
legacyVTKReader8Display = Show(legacyVTKReader8, renderView1)
# trace defaults for the display properties.
legacyVTKReader8Display.Representation = 'Surface With Edges'
legacyVTKReader8Display.ColorArrayName = ['POINTS', 'RegionId']
legacyVTKReader8Display.LookupTable = regionIdLUT
legacyVTKReader8Display.EdgeColor = [0.0, 0.0, 0.0]
legacyVTKReader8Display.OSPRayScaleArray = 'Normals'
legacyVTKReader8Display.OSPRayScaleFunction = 'PiecewiseFunction'
legacyVTKReader8Display.SelectOrientationVectors = 'None'
legacyVTKReader8Display.ScaleFactor = 3.95
legacyVTKReader8Display.SelectScaleArray = 'None'
legacyVTKReader8Display.GlyphType = 'Arrow'
legacyVTKReader8Display.GlyphTableIndexArray = 'None'
legacyVTKReader8Display.DataAxesGrid = 'GridAxesRepresentation'
legacyVTKReader8Display.PolarAxes = 'PolarAxesRepresentation'
legacyVTKReader8Display.GaussianRadius = 1.975
legacyVTKReader8Display.SetScaleArray = ['POINTS', '']
legacyVTKReader8Display.ScaleTransferFunction = 'PiecewiseFunction'
legacyVTKReader8Display.OpacityArray = ['POINTS', '']
legacyVTKReader8Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
legacyVTKReader8Display.SetScalarBarVisibility(renderView1, True)

# show data from transform5
transform5Display = Show(transform5, renderView1)
# trace defaults for the display properties.
transform5Display.Representation = 'Surface With Edges'
transform5Display.ColorArrayName = ['POINTS', 'RegionId']
transform5Display.LookupTable = regionIdLUT
transform5Display.EdgeColor = [0.0, 0.0, 0.0]
transform5Display.OSPRayScaleArray = 'Normals'
transform5Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform5Display.SelectOrientationVectors = 'None'
transform5Display.ScaleFactor = 0.47500000000000003
transform5Display.SelectScaleArray = 'None'
transform5Display.GlyphType = 'Arrow'
transform5Display.GlyphTableIndexArray = 'None'
transform5Display.DataAxesGrid = 'GridAxesRepresentation'
transform5Display.PolarAxes = 'PolarAxesRepresentation'
transform5Display.GaussianRadius = 0.23750000000000002
transform5Display.SetScaleArray = ['POINTS', '']
transform5Display.ScaleTransferFunction = 'PiecewiseFunction'
transform5Display.OpacityArray = ['POINTS', '']
transform5Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform5Display.SetScalarBarVisibility(renderView1, True)

# show data from transform6
transform6Display = Show(transform6, renderView1)
# trace defaults for the display properties.
transform6Display.Representation = 'Surface With Edges'
transform6Display.ColorArrayName = ['POINTS', 'RegionId']
transform6Display.LookupTable = regionIdLUT
transform6Display.EdgeColor = [0.0, 0.0, 0.0]
transform6Display.OSPRayScaleArray = 'Normals'
transform6Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform6Display.SelectOrientationVectors = 'None'
transform6Display.ScaleFactor = 4.383330154418945
transform6Display.SelectScaleArray = 'None'
transform6Display.GlyphType = 'Arrow'
transform6Display.GlyphTableIndexArray = 'None'
transform6Display.DataAxesGrid = 'GridAxesRepresentation'
transform6Display.PolarAxes = 'PolarAxesRepresentation'
transform6Display.GaussianRadius = 2.1916650772094726
transform6Display.SetScaleArray = ['POINTS', '']
transform6Display.ScaleTransferFunction = 'PiecewiseFunction'
transform6Display.OpacityArray = ['POINTS', '']
transform6Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform6Display.SetScalarBarVisibility(renderView1, True)

# show data from transform1
transform1Display = Show(transform1, renderView1)
# trace defaults for the display properties.
transform1Display.Representation = 'Surface With Edges'
transform1Display.ColorArrayName = ['POINTS', 'RegionId']
transform1Display.LookupTable = regionIdLUT
transform1Display.EdgeColor = [0.0, 0.0, 0.0]
transform1Display.OSPRayScaleArray = 'Normals'
transform1Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform1Display.SelectOrientationVectors = 'None'
transform1Display.ScaleFactor = 1.2875
transform1Display.SelectScaleArray = 'None'
transform1Display.GlyphType = 'Arrow'
transform1Display.GlyphTableIndexArray = 'None'
transform1Display.DataAxesGrid = 'GridAxesRepresentation'
transform1Display.PolarAxes = 'PolarAxesRepresentation'
transform1Display.GaussianRadius = 0.64375
transform1Display.SetScaleArray = ['POINTS', '']
transform1Display.ScaleTransferFunction = 'PiecewiseFunction'
transform1Display.OpacityArray = ['POINTS', '']
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
transform3Display.ScaleFactor = 0.9124997000000001
transform3Display.SelectScaleArray = 'RegionId'
transform3Display.GlyphType = 'Arrow'
transform3Display.GlyphTableIndexArray = 'RegionId'
transform3Display.DataAxesGrid = 'GridAxesRepresentation'
transform3Display.PolarAxes = 'PolarAxesRepresentation'
transform3Display.GaussianRadius = 0.45624985000000007
transform3Display.SetScaleArray = ['POINTS', 'RegionId']
transform3Display.ScaleTransferFunction = 'PiecewiseFunction'
transform3Display.OpacityArray = ['POINTS', 'RegionId']
transform3Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform3Display.SetScalarBarVisibility(renderView1, True)

# show data from transform7
transform7Display = Show(transform7, renderView1)
# trace defaults for the display properties.
transform7Display.Representation = 'Surface With Edges'
transform7Display.ColorArrayName = ['POINTS', 'RegionId']
transform7Display.LookupTable = regionIdLUT
transform7Display.EdgeColor = [0.0, 0.0, 0.0]
transform7Display.OSPRayScaleArray = 'RegionId'
transform7Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform7Display.SelectOrientationVectors = 'None'
transform7Display.ScaleFactor = 5.1141066
transform7Display.SelectScaleArray = 'RegionId'
transform7Display.GlyphType = 'Arrow'
transform7Display.GlyphTableIndexArray = 'RegionId'
transform7Display.DataAxesGrid = 'GridAxesRepresentation'
transform7Display.PolarAxes = 'PolarAxesRepresentation'
transform7Display.GaussianRadius = 2.5570533
transform7Display.SetScaleArray = ['POINTS', 'RegionId']
transform7Display.ScaleTransferFunction = 'PiecewiseFunction'
transform7Display.OpacityArray = ['POINTS', 'RegionId']
transform7Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform7Display.SetScalarBarVisibility(renderView1, True)

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
transform4Display.ScaleFactor = 1.8833566000000002
transform4Display.SelectScaleArray = 'RegionId'
transform4Display.GlyphType = 'Arrow'
transform4Display.GlyphTableIndexArray = 'RegionId'
transform4Display.DataAxesGrid = 'GridAxesRepresentation'
transform4Display.PolarAxes = 'PolarAxesRepresentation'
transform4Display.GaussianRadius = 0.9416783000000001
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
transform2Display.ScaleFactor = 0.2624997
transform2Display.SelectScaleArray = 'RegionId'
transform2Display.GlyphType = 'Arrow'
transform2Display.GlyphTableIndexArray = 'RegionId'
transform2Display.DataAxesGrid = 'GridAxesRepresentation'
transform2Display.PolarAxes = 'PolarAxesRepresentation'
transform2Display.GaussianRadius = 0.13124985
transform2Display.SetScaleArray = ['POINTS', 'RegionId']
transform2Display.ScaleTransferFunction = 'PiecewiseFunction'
transform2Display.OpacityArray = ['POINTS', 'RegionId']
transform2Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform2Display.SetScalarBarVisibility(renderView1, True)

# show data from transform8
transform8Display = Show(transform8, renderView1)
# trace defaults for the display properties.
transform8Display.Representation = 'Surface With Edges'
transform8Display.ColorArrayName = ['POINTS', 'RegionId']
transform8Display.LookupTable = regionIdLUT
transform8Display.EdgeColor = [0.0, 0.0, 0.0]
transform8Display.OSPRayScaleArray = 'RegionId'
transform8Display.OSPRayScaleFunction = 'PiecewiseFunction'
transform8Display.SelectOrientationVectors = 'None'
transform8Display.ScaleFactor = 1.2999967000000001
transform8Display.SelectScaleArray = 'RegionId'
transform8Display.GlyphType = 'Arrow'
transform8Display.GlyphTableIndexArray = 'RegionId'
transform8Display.DataAxesGrid = 'GridAxesRepresentation'
transform8Display.PolarAxes = 'PolarAxesRepresentation'
transform8Display.GaussianRadius = 0.6499983500000001
transform8Display.SetScaleArray = ['POINTS', 'RegionId']
transform8Display.ScaleTransferFunction = 'PiecewiseFunction'
transform8Display.OpacityArray = ['POINTS', 'RegionId']
transform8Display.OpacityTransferFunction = 'PiecewiseFunction'

# show color legend
transform8Display.SetScalarBarVisibility(renderView1, True)

# setup the color legend parameters for each legend in this view

# get color legend/bar for regionIdLUT in view renderView1
regionIdLUTColorBar = GetScalarBar(regionIdLUT, renderView1)
regionIdLUTColorBar.Title = 'RegionId'
regionIdLUTColorBar.ComponentTitle = ''

# ----------------------------------------------------------------
# finally, restore active source
SetActiveSource(transform8)
# ----------------------------------------------------------------
