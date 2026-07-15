import picocli.CommandLine;

import javax.swing.SwingUtilities;

public class Launcher {
    public static void main(String[] args) {
        if (args.length == 0) {
            SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                    new DunesGui().setVisible(true);
                }
            });
        } else {
            CommandLine.run(new Main(), args);
        }
    }
}
