// QuPath 0.7 GUI startup script. Apply once per image open, after QuPath restores its display.
import javafx.application.Platform
import javafx.collections.ListChangeListener
import qupath.lib.display.DirectServerChannelInfo
import qupath.lib.gui.QuPathGUI
import qupath.lib.gui.viewer.QuPathViewer
import qupath.lib.gui.viewer.QuPathViewerListener
import qupath.lib.images.ImageData
import qupath.lib.objects.PathObject

import java.awt.Shape
import java.awt.image.BufferedImage

class ReviewChannelListener implements QuPathViewerListener {
    @Override
    void imageDataChanged(QuPathViewer viewer, ImageData<BufferedImage> previous, ImageData<BufferedImage> current) {
        if (current == null || current.isBrightfield() || current.getServer().isRGB())
            return

        def display = viewer.getImageDisplay()
        def channels = display.availableChannels().findAll { it instanceof DirectServerChannelInfo }
        if (channels.size() < 2 || !display.selectedChannels().containsAll(channels))
            return

        // Multiplex panels may contain a nuclear stain for every cycle and blank exposures.
        def names = channels.collectEntries { [(it): it.getOriginalChannelName().toLowerCase(Locale.ROOT)] }
        def dapi = channels.find { names[it] == 'dapi' } ?: channels.find { names[it].contains('dapi') }
        def nuclear = dapi ?: channels.find { names[it].contains('hoechst') }
        def markers = channels.findAll {
            !(names[it] =~ /dapi|hoechst|blank|autofluorescen|(?:^|[^a-z0-9])(?:af|background)(?:$|[^a-z0-9])/)
        }
        Collections.shuffle(markers)
        def keep = nuclear == null ? [] : [nuclear]
        int markerCount = Math.min(3, channels.size() - 1 - keep.size())
        keep.addAll(markers.take(markerCount))
        if (keep.isEmpty())
            keep.add(channels.first())
        channels.findAll { !keep.contains(it) }.each { display.setChannelSelected(it, false) }
    }

    @Override void visibleRegionChanged(QuPathViewer viewer, Shape region) {}
    @Override void selectedObjectChanged(QuPathViewer viewer, PathObject object) {}
    @Override void viewerClosed(QuPathViewer viewer) {}
}

def install = {
    def gui = QuPathGUI.getInstance()
    def properties = gui.getStage().getProperties()
    def key = 'isletExplorer.defaultReviewChannels'
    if (properties.containsKey(key))
        return
    def listener = new ReviewChannelListener()
    def attach = { viewer ->
        viewer.addViewerListener(listener)
        listener.imageDataChanged(viewer, null, viewer.getImageData())
    }
    def viewers = gui.getAllViewers()
    viewers.each(attach)
    viewers.addListener({ change ->
        while (change.next()) {
            change.removed.each { it.removeViewerListener(listener) }
            change.addedSubList.each(attach)
        }
    } as ListChangeListener)
    properties.put(key, listener)
}
if (Platform.isFxApplicationThread())
    install()
else
    Platform.runLater(install as Runnable)
