// Download this file, open it in Fiji's Script Editor and press Run.
// The shared installer below handles downloads, verification, backup and replacement.
if (binding.hasVariable("CALCIUM_INSTALL_TEST")) return CalciumInstaller
def window = new javax.swing.JFrame("Install " + CalciumRelease.NAME)
def label = new javax.swing.JLabel("Preparing installation...")
def bar = new javax.swing.JProgressBar(0,100); bar.setStringPainted(true)
def close = new javax.swing.JButton("Close")
def panel = new javax.swing.JPanel(); panel.setLayout(new javax.swing.BoxLayout(panel,javax.swing.BoxLayout.Y_AXIS))
panel.setBorder(javax.swing.BorderFactory.createEmptyBorder(12,12,12,12))
label.setPreferredSize(new java.awt.Dimension(460,140)); label.setMaximumSize(new java.awt.Dimension(460,140))
panel.add(label); panel.add(bar); panel.add(close)
close.addActionListener({ e -> window.dispose() } as java.awt.event.ActionListener)
def show = {
    window.setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE)
    window.add(panel); window.pack(); window.setLocationRelativeTo(null); window.setVisible(true)
}
if (javax.swing.SwingUtilities.isEventDispatchThread()) show.call()
else javax.swing.SwingUtilities.invokeAndWait(show as Runnable)
Thread.start("Install calcium plugin") {
    try {
        def release = CalciumRelease.latest()
        String directory = ij.IJ.getDirectory("imagej")
        def installed = new CalciumInstaller(applicationDirectory: directory ? new File(directory) : null).install(release, { int percent, String message ->
            javax.swing.SwingUtilities.invokeLater { if(window.isDisplayable()) { bar.setValue(percent); bar.setString(message) } }
        })
        javax.swing.SwingUtilities.invokeLater {
            if (window.isDisplayable()) {
                bar.setValue(100); bar.setString("Installed " + installed.version)
                label.setText("<html><body style='width:340px'><b>Installed " + installed.version + ".</b><br>Save your work and restart Fiji.<br>Open Plugins &gt; " + CalciumInstaller.escape(CalciumRelease.command(installed.target.path)) + ".<br>Future releases can be installed with the plugin's Install update button." +
                    (installed.backup ? "<br>Your previous file is in calcium-plugin-backups inside Fiji's folder." : "") + "</body></html>")
            }
        }
    } catch (Exception problem) {
        javax.swing.SwingUtilities.invokeLater {
            if (window.isDisplayable()) {
                bar.setString("Installation did not finish")
                label.setText("<html><body style='width:340px'>Could not install: " + CalciumInstaller.escape(problem.message ?: problem.toString()) + "<br>Close this window and run the installer again to retry.</body></html>")
            }
        }
    }
}
