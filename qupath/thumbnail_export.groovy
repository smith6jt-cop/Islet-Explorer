// QuPath script: Export islet thumbnails (RGB) with overlays & manifest
// - Exports regions per 'Islet' using the matching 'IsletExpanded' (20 µm) when present
// - Draws outlines for classes: Islet, IsletExpanded, Nerve, Capillary, Lymphatic
// - Uses current viewer channel visibility/colors (turn on DAPI, INS, GCG, SST, CD3e, CD68, CD163, PDPN, CD31, B3TUBB, LGALS3 before running)
// - Writes a manifest: donor_id,islet_id,filename

import qupath.lib.gui.scripting.QPEx
import qupath.lib.objects.PathObject
import qupath.lib.regions.RegionRequest
// Note: We avoid RenderedImageServer to maximize compatibility across QuPath versions.
import java.awt.BasicStroke
import java.awt.Color
import java.awt.geom.AffineTransform
import java.awt.image.BufferedImage
import javax.imageio.ImageIO
import javax.swing.JFileChooser

// -------- configure --------
// Fixed output directory per request
def OUTPUT_DIR = new File('/home/smith6jt/panc_CODEX/streamlit_app/qupath')
OUTPUT_DIR.mkdirs()

def DONOR_ID = null // set explicitly or leave null to use image name
double DOWNSAMPLE = 4.0 // Higher -> smaller images (4–12 typical); 4 gives higher resolution

// -------- helpers --------
def getIsletId = { PathObject obj ->
    def ml = obj.getMeasurementList()
    def id = null
    for (key in ['islet_id', 'Islet_ID', 'isletId', 'ID', 'Name']) {
        if (ml != null && ml.containsKey(key)) {
            id = ml.getMeasurementValue(key)
            break
        }
    }
    if (id == null) id = obj.getName()
    return (id == null) ? null : id.toString()
}

def classColor = { String cls ->
    switch (cls) {
        case 'Islet': return new Color(0, 255, 0)        // green
        case 'IsletExpanded': return new Color(0, 150, 0) // dark green
        case 'Nerve': return Color.YELLOW
        case 'Capillary': return new Color(255, 80, 80)  // reddish
        case 'Lymphatic': return Color.CYAN
        default: return Color.WHITE
    }
}

// -------- core --------
def imageData = QPEx.getCurrentImageData()
def server = imageData.getServer()
def imageName = server.getMetadata().getName()
def donorId = DONOR_ID ?: imageName

def allAnn = QPEx.getAnnotationObjects()

def islets = allAnn.findAll { it.getPathClass() != null && it.getPathClass().toString().equalsIgnoreCase('Islet') && it.getROI() != null }
def expanded = allAnn.findAll { it.getPathClass() != null && it.getPathClass().toString().equalsIgnoreCase('IsletExpanded') && it.getROI() != null }

def others = allAnn.findAll { it.getROI() != null && it.getPathClass() != null && ['Nerve','Capillary','Lymphatic'].contains(it.getPathClass().toString()) }

// Prepare output
def outDir = (OUTPUT_DIR instanceof File) ? OUTPUT_DIR : new File(OUTPUT_DIR as String)
outDir.mkdirs()
def manifestFile = new File(outDir, 'manifest.csv')
def manifestHeader = !manifestFile.exists()
def writer = new java.io.FileWriter(manifestFile, true)
if (manifestHeader) writer.write('donor_id,islet_id,filename\n')

// Note: We render regions directly from the image server for maximum compatibility.
// This won't apply viewer channel blending, but overlays will still be drawn.

int exported = 0
for (islet in islets) {
    def roi = islet.getROI()
    if (roi == null) continue

    // Find matching expanded object by overlap
    def exp = expanded.find { it.getROI() != null && it.getROI().getGeometry().intersects(roi.getGeometry()) }
    def exportRoi = (exp != null) ? exp.getROI() : roi

    def isletId = getIsletId(islet) ?: (exported + 1).toString()
    def safeIsletId = isletId.replaceAll('[^A-Za-z0-9._-]', '_')
    def safeDonor = donorId.replaceAll('[^A-Za-z0-9._-]', '_')
    def filename = String.format('%s_%s.png', safeDonor, safeIsletId)
    def outFile = new File(outDir, filename)

    try {
        def req = RegionRequest.createInstance(server.getPath(), DOWNSAMPLE, exportRoi)
        def imgSrc = server.readRegion(req)

        // Copy to ARGB to ensure a compatible color model for drawing
        def out = new BufferedImage(imgSrc.getWidth(), imgSrc.getHeight(), BufferedImage.TYPE_INT_ARGB)
        def g2d = out.createGraphics()
        g2d.drawImage(imgSrc, 0, 0, null)
        g2d.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
        def at = new AffineTransform()
        // Translate by top-left of the requested region in full-resolution pixels,
        // then scale down to match the requested downsample.
        at.translate(-req.getX() as double, -req.getY() as double)
        at.scale(1.0 / DOWNSAMPLE, 1.0 / DOWNSAMPLE)
        g2d.transform(at)

        def drawObj = { PathObject obj, float strokeWidth, boolean dashed ->
            def shape = obj.getROI().getShape()
            def col = classColor(obj.getPathClass().toString())
            // Fill with translucent color to emphasize annotations
            def fill = new Color(col.getRed(), col.getGreen(), col.getBlue(), 60)
            g2d.setPaint(fill)
            g2d.fill(shape)
            g2d.setColor(col)
            def stroke = dashed ? new BasicStroke(strokeWidth, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND, 10f, [10f, 10f] as float[], 0f) : new BasicStroke(strokeWidth)
            g2d.setStroke(stroke)
            g2d.draw(shape)
        }

        // Draw Islet & IsletExpanded first
        drawObj(islet, 4f, false)
        if (exp != null) drawObj(exp, 3f, true)
        // Draw other objects within the export ROI bounds
        for (o in others) {
            try {
                if (o.getROI().getBoundsX() > exportRoi.getBoundsX() + exportRoi.getBoundsWidth()) continue
                if (o.getROI().getBoundsY() > exportRoi.getBoundsY() + exportRoi.getBoundsHeight()) continue
                if (o.getROI().getBoundsX() + o.getROI().getBoundsWidth() < exportRoi.getBoundsX()) continue
                if (o.getROI().getBoundsY() + o.getROI().getBoundsHeight() < exportRoi.getBoundsY()) continue
                drawObj(o, 2f, false)
            } catch (Exception ex) {
                // ignore
            }
        }
        g2d.dispose()

        ImageIO.write(out, 'PNG', outFile)
        writer.write(String.format('%s,%s,%s\n', safeDonor, safeIsletId, filename))
        exported++
    } catch (Exception e) {
        print("Failed to export thumbnail for islet ${isletId}: ${e}")
    }
}

writer.close()
print("Exported " + exported + " thumbnails to " + outDir.getAbsolutePath())
