import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO

// The original 6090 acquisition is pixel-matched to its source project entry.
// Recovered 6539 uses that panel's channel colors and its own recorded contrast.
def aliases = ['DAPI-01':'DAPI', COL4A1:'ColIV', ATP1A1:'NaKATPase', CTNNB1:'BCatenin', ACTB:'BActin', KRT19:'CK19']
def names = getCurrentServer().getMetadata().getChannels()*.getName()
def normalized = names.collect { name ->
    String clean = name.replaceAll('^reg\\d+_cyc\\d+_ch\\d+_', '').replaceAll('(?i)\\.(tif|tiff)\\s*$', '').trim()
    aliases.get(clean, clean)
}
if (normalized.toSet().size() != normalized.size()) throw new IllegalStateException('Duplicate marker names')
String caseId = System.getenv('ISLET_CASE_ID')
if (!(caseId in ['6090','6539'])) throw new IllegalStateException('Unexpected recovered case')
def propertyName = 'qupath.lib.display.ImageDisplay'
def source = PathIO.readImageData(new File('/data/QuPath_Projects/QuPath_project_IsletExplorer/data/39/data.qpdata'))
def display = GsonTools.getInstance().fromJson(source.getProperty(propertyName) as String, List.class)
if (display.size() != normalized.size()) throw new IllegalStateException('Display channel count mismatch')
def ranges = caseId == '6539' ? GsonTools.getInstance().fromJson(
    new File(System.getenv('ISLET_EXPORT_DIR'), '../../display_ranges_6539.json').text, Map.class).channels : null
display.eachWithIndex { channel, i ->
    channel.name = "${normalized[i]} (C${i + 1})".toString()
    if (ranges != null) {
        channel.minDisplay = ranges[i].minDisplay
        channel.maxDisplay = ranges[i].maxDisplay
    }
}
setChannelNames(*normalized as String[])
setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.3774, 0.3774)
getCurrentImageData().setProperty(propertyName, GsonTools.getInstance().toJson(display))
println 'Recovered acquisition channels: ' + normalized
println 'Preserved source-project annotations separately; new image hierarchy has ' + getAnnotationObjects().size() + ' annotations'
