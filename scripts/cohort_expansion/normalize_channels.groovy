import qupath.lib.io.GsonTools
import qupath.lib.io.PathIO

def aliases = [DAPI_01:'DAPI', 'DAPI-01':'DAPI', COL4A1:'ColIV', ATP1A1:'NaKATPase', CTNNB1:'BCatenin', ACTB:'BActin', KRT19:'CK19']
def originalNames = getCurrentServer().getMetadata().getChannels()*.getName()
def normalized = originalNames.collect { name ->
    String clean = name.replaceAll('^reg\\d+_cyc\\d+_ch\\d+_', '').replaceAll('\\s*\\(C\\d+\\)\\s*$', '').replaceAll('(?i)\\.(tif|tiff)\\s*$', '').trim()
    aliases.get(clean, clean)
}
if (normalized.contains('KRT8')) {
    if (System.getenv('ISLET_KRT8_PANEL_VERIFIED') != '1')
        throw new IllegalStateException('KRT8 acquisition equivalence has not been verified')
    normalized = normalized.collect { it == 'KRT8' ? 'KRT8-18' : it }
}
if (normalized.toSet().size() != normalized.size()) throw new IllegalStateException('Duplicate canonical marker names')
def propertyName = 'qupath.lib.display.ImageDisplay'
def sourceEntry = [6232:43, 6558:44, 6520:40, 78:42][Integer.parseInt(System.getenv('ISLET_CASE_ID'))]
def sourceData = PathIO.readImageData(new File('/data/QuPath_Projects/QuPath_project_IsletExplorer/data/' + sourceEntry + '/data.qpdata'))
def saved = sourceData.getProperty(propertyName)
if (saved) {
    def display = GsonTools.getInstance().fromJson(saved as String, List.class)
    display.eachWithIndex { channel, i -> channel.name = "${normalized[i]} (C${i + 1})".toString() }
    getCurrentImageData().setProperty(propertyName, GsonTools.getInstance().toJson(display))
}
setChannelNames(*normalized as String[])
setImageType('FLUORESCENCE')
setPixelSizeMicrons(0.3774, 0.3774)
println 'Canonical channels: ' + normalized
println 'Source annotations preserved until tissue/islet stage: ' + getAnnotationObjects().size()
