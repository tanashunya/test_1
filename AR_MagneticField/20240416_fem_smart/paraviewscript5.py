import config
# trace generated using paraview version 5.10.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
import os

def makevecpng_DL(current_dir, filename):
    # # trace generated using paraview version 5.10.0-RC1
    #import paraview
    #paraview.compatibility.major = 5
    #paraview.compatibility.minor = 10

    #### import the simple module from the paraview
    # from paraview.simple import *
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Legacy VTK Reader'
    dL_vectorvtk = LegacyVTKReader(registrationName='DL_vector.vtk', FileNames=[current_dir+'/DL_vector.vtk'])

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
    dL_vectorvtkDisplay.ScaleFactor = 64.0
    dL_vectorvtkDisplay.SelectScaleArray = 'None'
    dL_vectorvtkDisplay.GlyphType = 'Arrow'
    dL_vectorvtkDisplay.GlyphTableIndexArray = 'None'
    dL_vectorvtkDisplay.GaussianRadius = 3.2
    dL_vectorvtkDisplay.SetScaleArray = [None, '']
    dL_vectorvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    dL_vectorvtkDisplay.OpacityArray = [None, '']
    dL_vectorvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    dL_vectorvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
    dL_vectorvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
    dL_vectorvtkDisplay.ScalarOpacityUnitDistance = 0
    dL_vectorvtkDisplay.ScalarOpacityFunction = bPWF
    dL_vectorvtkDisplay.OpacityArrayName = ['CELLS', 'B']
    dL_vectorvtkDisplay.SliceFunction = 'Plane'

    # init the 'Plane' selected for 'SliceFunction'
    dL_vectorvtkDisplay.SliceFunction.Origin = [128.0, 128.0, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [128.0, 128.0, -10000.0]
    renderView1.CameraFocalPoint = [128.0, 128.0, 0.0]

    # show color bar/color legend
    dL_vectorvtkDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    bLUT.ApplyPreset('Rainbow Uniform', True)

    # create a new 'Glyph'
    glyph1 = Glyph(registrationName='Glyph1', Input=dL_vectorvtk,
        GlyphType='Arrow')
    glyph1.OrientationArray = ['CELLS', 'B']
    glyph1.ScaleArray = ['POINTS', 'No scale array']
    glyph1.ScaleFactor = 20
    glyph1.GlyphTransform = 'Transform2'

    # set active source
    SetActiveSource(dL_vectorvtk)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Show3DWidgets(proxy=dL_vectorvtkDisplay.SliceFunction)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Show3DWidgets(proxy=dL_vectorvtkDisplay)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Hide3DWidgets(proxy=dL_vectorvtkDisplay.SliceFunction)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Hide3DWidgets(proxy=dL_vectorvtkDisplay)

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
    glyph1Display.ScaleFactor = 63.29572865962982
    glyph1Display.SelectScaleArray = 'B'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'B'
    glyph1Display.GaussianRadius = 3.164786432981491
    glyph1Display.SetScaleArray = ['POINTS', 'B']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'B']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'
    # Properties modified on glyph1
    glyph1.GlyphMode = 'All Points'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    # glyph1Display.ScaleTransferFunction.Points = [-0.29878848791122437, 0.0, 0.5, 0.0, 0.6960177421569824, 1.0, 0.5, 0.0]
    glyph1Display.ScaleTransferFunction.Points = [-0.2640950083732605, 0.0, 0.5, 0.0, 0.39233338832855225, 1.0, 0.5, 0.0]
    

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    # glyph1Display.OpacityTransferFunction.Points = [-0.29878848791122437, 0.0, 0.5, 0.0, 0.6960177421569824, 1.0, 0.5, 0.0]
    glyph1Display.OpacityTransferFunction.Points = [-0.2640950083732605, 0.0, 0.5, 0.0, 0.39233338832855225, 1.0, 0.5, 0.0]

    # set scalar coloring
    ColorBy(glyph1Display, ('POINTS', 'B', 'Magnitude'))

    # rescale color and/or opacity maps used to include current data range
    glyph1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, config.COLOR_BAR)

    # get color transfer function/color map for 'cell_data'
    cell_dataLUT = GetColorTransferFunction('B')

    # Properties modified on glyph1
    glyph1.ScaleFactor = 23.68


    # Rescale transfer function
    bLUT.RescaleTransferFunction(0.0, 1.3)

    # Rescale transfer function
    bPWF.RescaleTransferFunction(0.0, 1.3)

    # reset view to fit data
    renderView1.ResetCamera(False)

    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0

    # hide data in view
    Hide(glyph1, renderView1)

    # reset view to fit data
    renderView1.ResetCamera(False)

    # Properties modified on glyph1
    glyph1.GlyphMode = 'All Points'

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # update the view to ensure updated data information
    renderView1.Update()




    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    cell_dataLUT.ApplyPreset('Rainbow Uniform', True)

    # get color legend/bar for cell_dataLUT in view renderView1
    cell_dataLUTColorBar = GetScalarBar(cell_dataLUT, renderView1)

    # change scalar bar placement
    cell_dataLUTColorBar.Position = [0.8265625, 0.025]

    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0

    # Properties modified on cell_dataLUTColorBar
    cell_dataLUTColorBar.Title = 'B [T]'
    cell_dataLUTColorBar.ComponentTitle = ''
    cell_dataLUTColorBar.TitleFontFamily = 'Times'
    cell_dataLUTColorBar.LabelFontFamily = 'Times'

    # hide data in view
    Hide(dL_vectorvtk, renderView1)

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(256, 256)

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [128.0, 128.0, -1]
    renderView1.CameraFocalPoint = [128.0, 128.0, 0.0]
    renderView1.CameraParallelScale = 225.7895720215109 * config.HEIGHT / 480
    # renderView1.Update()

    # Properties modified on glyph1
    glyph1.GlyphMode = 'All Points'

    # save screenshot
    SaveScreenshot(current_dir+'/paraview_img_DL.png', renderView1, ImageResolution=[256, 256],
        TransparentBackground=1)

    # # set active source
    # SetActiveSource(dL_vectorvtk)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Show3DWidgets(proxy=dL_vectorvtkDisplay.SliceFunction)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Show3DWidgets(proxy=dL_vectorvtkDisplay)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Hide3DWidgets(proxy=dL_vectorvtkDisplay.SliceFunction)

    # # toggle 3D widget visibility (only when running from the GUI)
    # Hide3DWidgets(proxy=dL_vectorvtkDisplay)

    # # hide data in view
    # Hide(glyph1, renderView1)

    # show data in view
    # dL_vectorvtkDisplay = Show(dL_vectorvtk, renderView1, 'UniformGridRepresentation')

    # # show color bar/color legend
    # dL_vectorvtkDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy glyph1
    Delete(glyph1)
    del glyph1

    # destroy dL_vectorvtk
    Delete(dL_vectorvtk)
    del dL_vectorvtk

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(256, 256)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [128.0, 128.0, -1]
    renderView1.CameraFocalPoint = [128.0, 128.0, 0.0]
    renderView1.CameraParallelScale = 400

    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).









def makevecpng(current_dir):
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Legacy VTK Reader'
    # print(current_directory)
    vectorsvtk = LegacyVTKReader(registrationName='vectors.vtk', FileNames=[current_dir+'/vectors.vtk'])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    vectorsvtkDisplay = Show(vectorsvtk, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'B'
    bLUT = GetColorTransferFunction('B')

    # get opacity transfer function/opacity map for 'B'
    bPWF = GetOpacityTransferFunction('B')

    # trace defaults for the display properties.
    vectorsvtkDisplay.Representation = 'Surface'
    vectorsvtkDisplay.ColorArrayName = [None, '']
    vectorsvtkDisplay.SelectTCoordArray = 'None'
    vectorsvtkDisplay.SelectNormalArray = 'None'
    vectorsvtkDisplay.SelectTangentArray = 'None'
    vectorsvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    vectorsvtkDisplay.SelectOrientationVectors = 'cell_data'
    vectorsvtkDisplay.ScaleFactor = 64.0
    vectorsvtkDisplay.SelectScaleArray = 'None'
    vectorsvtkDisplay.GlyphType = 'Arrow'
    vectorsvtkDisplay.GlyphTableIndexArray = 'None'
    vectorsvtkDisplay.GaussianRadius = 3.2
    vectorsvtkDisplay.SetScaleArray = [None, '']
    vectorsvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    vectorsvtkDisplay.OpacityArray = [None, '']
    vectorsvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    vectorsvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
    vectorsvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
    vectorsvtkDisplay.ScalarOpacityUnitDistance = 0
    vectorsvtkDisplay.OpacityArrayName = ['CELLS', 'cell_data']

    # reset view to fit data
    renderView1.ResetCamera(False)

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 10000.0]
    renderView1.CameraFocalPoint = [config.WIDTH/2.0, config.HEIGHT/2.0, 0.0]

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Glyph'
    glyph1 = Glyph(registrationName='Glyph1', Input=vectorsvtk,
        GlyphType='Arrow')
    glyph1.OrientationArray = ['CELLS', 'cell_data']
    glyph1.ScaleArray = ['POINTS', 'No scale array']
    glyph1.ScaleFactor = 64.0
    glyph1.GlyphTransform = 'Transform2'

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
    glyph1Display.OSPRayScaleArray = 'cell_data'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'cell_data'
    glyph1Display.ScaleFactor = 63.29572865962982
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = 3.164786432981491
    glyph1Display.SetScaleArray = ['POINTS', 'cell_data']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'cell_data']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    glyph1Display.ScaleTransferFunction.Points = [-0.2640950083732605, 0.0, 0.5, 0.0, 0.39233338832855225, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    glyph1Display.OpacityTransferFunction.Points = [-0.2640950083732605, 0.0, 0.5, 0.0, 0.39233338832855225, 1.0, 0.5, 0.0]

    # hide data in view
    Hide(vectorsvtk, renderView1)

    # set scalar coloring
    ColorBy(glyph1Display, ('POINTS', 'cell_data', 'Magnitude'))

    # rescale color and/or opacity maps used to include current data range
    glyph1Display.RescaleTransferFunctionToDataRange(True, False)

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, config.COLOR_BAR)

    # get color transfer function/color map for 'cell_data'
    cell_dataLUT = GetColorTransferFunction('cell_data')

    # get opacity transfer function/opacity map for 'cell_data'
    # cell_dataPWF = GetOpacityTransferFunction('cell_data')

    # Properties modified on glyph1
    glyph1.ScaleFactor = 23.68

    # Rescale transfer function
    bLUT.RescaleTransferFunction(0.0, 1.3)

    # Rescale transfer function
    bPWF.RescaleTransferFunction(0.0, 1.3)

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # reset view to fit data
    renderView1.ResetCamera(False)

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '3D'

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, config.COLOR_BAR)



    # update the view to ensure updated data information
    renderView1.Update()

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    cell_dataLUT.ApplyPreset('Rainbow Uniform', True)

    # get color legend/bar for cell_dataLUT in view renderView1
    cell_dataLUTColorBar = GetScalarBar(cell_dataLUT, renderView1)

    # change scalar bar placement
    cell_dataLUTColorBar.Position = [0.8265625, 0.025]

    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0

    # Properties modified on cell_dataLUTColorBar
    cell_dataLUTColorBar.Title = 'B [T]'
    cell_dataLUTColorBar.ComponentTitle = ''
    cell_dataLUTColorBar.TitleFontFamily = 'Times'
    cell_dataLUTColorBar.LabelFontFamily = 'Times'

    # set active source
    SetActiveSource(vectorsvtk)

    # set active source
    # SetActiveSource(vectorsvtk)

    # show data in view
    vectorsvtkDisplay = Show(vectorsvtk, renderView1, 'UnstructuredGridRepresentation')

    # hide data in view
    Hide(glyph1, renderView1)
    

    # reset view to fit data
    renderView1.ResetCamera(False)

    # set active source
    SetActiveSource(glyph1)

    # show data in view
    glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)

    # Properties modified on glyph1
    glyph1.GlyphMode = 'All Points'

    # hide data in view
    Hide(vectorsvtk, renderView1)

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(config.WIDTH, config.HEIGHT)

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    # renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 892.2007692615181]
    renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 10000.0]
    # renderView1.CameraFocalPoint = [config.WIDTH/2.0, config.HEIGHT/2.0, -653.280552800991]
    renderView1.CameraFocalPoint = [config.WIDTH/2.0, config.HEIGHT/2.0, 0.0]
    renderView1.CameraParallelScale = 225.7895720215109 * config.HEIGHT / 480
    
    # if config.WIDTH==256 and config.HEIGHT==256:
    #     renderView1.CameraParallelScale = 106.01626497342944
    # else:
    #     renderView1.CameraParallelScale = 225.7895720215109 * config.WIDTH / 640

    glyph1Display.SetScalarBarVisibility(renderView1, config.COLOR_BAR)

    # save screenshot
    SaveScreenshot(current_dir+'/paraview_img.png', renderView1, ImageResolution=[config.WIDTH, config.HEIGHT],
        TransparentBackground=1)

    # set active source
    SetActiveSource(vectorsvtk)

    # hide data in view
    Hide(glyph1, renderView1)

    # show data in view
    vectorsvtkDisplay = Show(vectorsvtk, renderView1, 'UnstructuredGridRepresentation')

    # destroy glyph1
    Delete(glyph1)
    del glyph1

    # destroy vectorsvtk
    Delete(vectorsvtk)
    del vectorsvtk

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(config.WIDTH, config.HEIGHT)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 892.2007692615181]
    renderView1.CameraFocalPoint = [config.WIDTH/2.0, config.HEIGHT/2.0, -653.280552800991]
    renderView1.CameraParallelScale = 400.0

    #--------------------------------------------
    # uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).






















#######################################################################################################################################################################
def makeconpng(current_dir):
    # trace generated using paraview version 5.10.0-RC1
    #import paraview
    #paraview.compatibility.major = 5
    #paraview.compatibility.minor = 10

    #### import the simple module from the paraview
    # from paraview.simple import *
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'Legacy VTK Reader'
    contourvtk = LegacyVTKReader(registrationName='contour.vtk', FileNames=[current_dir+'/contour.vtk'])

    # set active source
    SetActiveSource(contourvtk)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # show data in view
    contourvtkDisplay = Show(contourvtk, renderView1, 'UnstructuredGridRepresentation')

    # get color transfer function/color map for 'A'
    aLUT = GetColorTransferFunction('A')

    # get opacity transfer function/opacity map for 'A'
    aPWF = GetOpacityTransferFunction('A')

    # trace defaults for the display properties.
    contourvtkDisplay.Representation = 'Surface'
    contourvtkDisplay.ColorArrayName = ['POINTS', 'A']
    contourvtkDisplay.LookupTable = aLUT
    contourvtkDisplay.SelectTCoordArray = 'None'
    contourvtkDisplay.SelectNormalArray = 'None'
    contourvtkDisplay.SelectTangentArray = 'None'
    contourvtkDisplay.OSPRayScaleArray = 'A'
    contourvtkDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    contourvtkDisplay.SelectOrientationVectors = 'cell_data'
    contourvtkDisplay.ScaleFactor = 64.0
    contourvtkDisplay.SelectScaleArray = 'A'
    contourvtkDisplay.GlyphType = 'Arrow'
    contourvtkDisplay.GlyphTableIndexArray = 'A'
    contourvtkDisplay.GaussianRadius = 3.2
    contourvtkDisplay.SetScaleArray = ['POINTS', 'A']
    contourvtkDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    contourvtkDisplay.OpacityArray = ['POINTS', 'A']
    contourvtkDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    contourvtkDisplay.DataAxesGrid = 'GridAxesRepresentation'
    contourvtkDisplay.PolarAxes = 'PolarAxesRepresentation'
    contourvtkDisplay.ScalarOpacityFunction = aPWF
    contourvtkDisplay.ScalarOpacityUnitDistance = 63.9904455201293
    contourvtkDisplay.OpacityArrayName = ['POINTS', 'A']

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contourvtkDisplay.ScaleTransferFunction.Points = [-22.99631690979004, 0.0, 0.5, 0.0, 27.742494583129883, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contourvtkDisplay.OpacityTransferFunction.Points = [-22.99631690979004, 0.0, 0.5, 0.0, 27.742494583129883, 1.0, 0.5, 0.0]

    # show color bar/color legend
    contourvtkDisplay.SetScalarBarVisibility(renderView1, True)

    #changing interaction mode based on data extents
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 10000.0]
    renderView1.CameraFocalPoint = [config.WIDTH/2.0, config.HEIGHT/2.0, 0.0]

    # reset view to fit data
    renderView1.ResetCamera(False)

    # show data in view
    contourvtkDisplay = Show(contourvtk, renderView1, 'UnstructuredGridRepresentation')

    # reset view to fit data
    renderView1.ResetCamera(False)

    #changing interaction mode based on data extents
    renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 10000.0]

    # show color bar/color legend
    contourvtkDisplay.SetScalarBarVisibility(renderView1, True)

    # update the view to ensure updated data information
    renderView1.Update()

    # create a new 'Contour'
    contour1 = Contour(registrationName='Contour1', Input=contourvtk)
    contour1.ContourBy = ['POINTS', 'A']
    contour1.Isosurfaces = [2.373088836669922]
    contour1.PointMergeMethod = 'Uniform Binning'

    # set active source
    SetActiveSource(contour1)

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = ['POINTS', 'A']
    contour1Display.LookupTable = aLUT
    contour1Display.SelectTCoordArray = 'None'
    contour1Display.SelectNormalArray = 'None'
    contour1Display.SelectTangentArray = 'None'
    contour1Display.OSPRayScaleArray = 'A'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'cell_data'
    contour1Display.ScaleFactor = 36.40548400878907
    contour1Display.SelectScaleArray = 'A'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'A'
    contour1Display.GaussianRadius = 1.8202742004394532
    contour1Display.SetScaleArray = ['POINTS', 'A']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'A']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.PolarAxes = 'PolarAxesRepresentation'

    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [2.373088836669922, 0.0, 0.5, 0.0, 2.373577117919922, 1.0, 0.5, 0.0]

    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [2.373088836669922, 0.0, 0.5, 0.0, 2.373577117919922, 1.0, 0.5, 0.0]

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)

    if config.DARK_CONTOUR:
        # turn off scalar coloring
        # ColorBy(contour1Display, None)

        # Hide the scalar bar for this color map if no visible data is colored by it.
        # HideScalarBarIfNotNeeded(aLUT, renderView1)

        # hide data in view
        # Hide(contourvtk, renderView1)

        # change solid color
        contour1Display.AmbientColor = [0.0, 0.0, 0.0]
        contour1Display.DiffuseColor = [0.0, 0.0, 0.0]

        # Properties modified on contour1
        # contour1.Isosurfaces = [3.9133596420288086, -15.166906356811523, -10.926847245958117, -6.68678813510471, -2.4467290242513027, 1.7933300866021042, 6.033389197455513, 10.273448308308918, 14.513507419162323, 18.753566530015732, 22.99362564086914]


    # Properties modified on contour1
    contour1.Isosurfaces = [2.373088836669922, -22.99631690979004, -20.32585314700478, -17.65538938421952, -14.984925621434261, -12.314461858649002, -9.643998095863743, -6.973534333078483, -4.303070570293226, -1.6326068075079654, 1.0378569552772952, 3.708320718062552, 6.378784480847813, 9.049248243633073, 11.719712006418334, 14.390175769203587, 17.060639531988848, 19.73110329477411, 22.40156705755937, 25.07203082034463, 27.742494583129883]

    # show data in view
    contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

    # hide data in view
    Hide(contourvtk, renderView1)

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)
    
    # update the view to ensure updated data information
    renderView1.Update()

    # Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
    aLUT.ApplyPreset('Rainbow Uniform', True)

    # Hide orientation axes
    renderView1.OrientationAxesVisibility = 0

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(config.WIDTH, config.HEIGHT)

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [config.WIDTH/2.0, config.HEIGHT/2.0, 10000.0]
    renderView1.CameraFocalPoint = [config.WIDTH/2.0, config.HEIGHT/2.0, 0.0]
    renderView1.CameraParallelScale = 225.7895720215109 * config.HEIGHT / 480

    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, config.COLOR_BAR)

    # save screenshot
    SaveScreenshot(current_dir+'/contour.png', renderView1, ImageResolution=[config.WIDTH, config.HEIGHT],
        TransparentBackground=1)

    # set active source
    SetActiveSource(contourvtk)

    # hide data in view
    Hide(contour1, renderView1)

    # show data in view
    contourvtkDisplay = Show(contourvtk, renderView1, 'UnstructuredGridRepresentation')

    # show color bar/color legend
    contourvtkDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy contour1
    Delete(contour1)
    del contour1

    # destroy contourvtk
    Delete(contourvtk)
    del contourvtk

    #================================================================
    # addendum: following script captures some of the application
    # state to faithfully reproduce the visualization during playback
    #================================================================

    #--------------------------------
    # saving layout sizes for layouts

    # layout/tab size in pixels
    layout1.SetSize(config.WIDTH, config.HEIGHT)

    #-----------------------------------
    # saving camera placements for views

    # current camera placement for renderView1
    renderView1.CameraPosition = [320.0, 240.0, 892.2007692615181]
    renderView1.CameraFocalPoint = [320.0, 240.0, -653.280552800991]
    renderView1.CameraParallelScale = 400.0