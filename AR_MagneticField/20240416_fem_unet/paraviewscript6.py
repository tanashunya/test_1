# trace generated using paraview version 5.10.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
dL_vectorvtk = LegacyVTKReader(registrationName='DL_vector.vtk', FileNames=['/home/syunya-linux/デスクトップ/MyWorkspace/20231120_write_vtk/DL_vector.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
dL_vectorvtkDisplay = Show(dL_vectorvtk, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'B'
bLUT = GetColorTransferFunction('B')

# get opacity transfer function/opacity map for 'B'
bPWF = GetOpacityTransferFunction('B')

# trace defaults for the display properties.
dL_vectorvtkDisplay.Representation = 'Slice'
dL_vectorvtkDisplay.ColorArrayName = ['CELLS', 'B']
dL_vectorvtkDisplay.LookupTable = bLUT
dL_vectorvtkDisplay.SelectTCoordArray = 'None'
dL_vectorvtkDisplay.SelectNormalArray = 'None'
dL_vectorvtkDisplay.SelectTangentArray = 'None'
dL_vectorvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
dL_vectorvtkDisplay.SelectOrientationVectors = 'B'
dL_vectorvtkDisplay.ScaleFactor = 25.6
dL_vectorvtkDisplay.SelectScaleArray = 'None'
dL_vectorvtkDisplay.GlyphType = 'Arrow'
dL_vectorvtkDisplay.GlyphTableIndexArray = 'None'
dL_vectorvtkDisplay.GaussianRadius = 1.28
dL_vectorvtkDisplay.SetScaleArray = [None, '']
dL_vectorvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
dL_vectorvtkDisplay.OpacityArray = [None, '']
dL_vectorvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
dL_vectorvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
dL_vectorvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
dL_vectorvtkDisplay.ScalarOpacityUnitDistance = 35.91878554589993
dL_vectorvtkDisplay.ScalarOpacityFunction = bPWF
dL_vectorvtkDisplay.OpacityArrayName = ['CELLS', 'B']
dL_vectorvtkDisplay.SliceFunction = 'Plane'

# init the 'Plane' selected for 'SliceFunction'
dL_vectorvtkDisplay.SliceFunction.Origin = [128.0, 128.0, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [128.0, 128.0, 10000.0]
renderView1.CameraFocalPoint = [128.0, 128.0, 0.0]

# show color bar/color legend
dL_vectorvtkDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=dL_vectorvtk,
    GlyphType='Arrow')
glyph1.OrientationArray = ['CELLS', 'B']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 25.6
glyph1.GlyphTransform = 'Transform2'

# reset view to fit data
renderView1.ResetCamera(False)

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# set active source
SetActiveSource(glyph1)

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = [None, '']
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'B'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'B'
glyph1Display.ScaleFactor = 29.627755546569826
glyph1Display.SelectScaleArray = 'B'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'B'
glyph1Display.GaussianRadius = 1.4813877773284911
glyph1Display.SetScaleArray = ['POINTS', 'B']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'B']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-0.2816946506500244, 0.0, 0.5, 0.0, 0.687757134437561, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-0.2816946506500244, 0.0, 0.5, 0.0, 0.687757134437561, 1.0, 0.5, 0.0]

# set scalar coloring
ColorBy(glyph1Display, ('POINTS', 'B', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
glyph1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
bLUT.ApplyPreset('Rainbow Uniform', True)

# hide data in view
Hide(dL_vectorvtk, renderView1)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(512, 512)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [128.0, 128.0, -699.4050067376327]
renderView1.CameraFocalPoint = [128.0, 128.0, 0.0]
renderView1.CameraParallelScale = 124.24226488854562

# save screenshot
SaveScreenshot('/home/syunya-linux/デスクトップ/MyWorkspace/20231120_write_vtk/paraview_img_DL.png', renderView1, ImageResolution=[256, 256],
    TransparentBackground=1, 
    # PNG options
    CompressionLevel='0')

# set active source
SetActiveSource(glyph1)

# set active source
SetActiveSource(None)

# set active view
SetActiveView(None)

# get active view
renderView1_1 = GetActiveViewOrCreate('RenderView')

# Create a new 'Render View'
renderView1_2 = CreateView('RenderView')
renderView1_2.AxesGrid = 'GridAxes3DActor'
renderView1_2.StereoType = 'Crystal Eyes'
renderView1_2.CameraFocalDisk = 1.0

# add view to a layout so it's visible in UI
AssignViewToLayout(view=renderView1_2, layout=None, hint=0)

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1_1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1_1.SetSize(512, 512)

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).