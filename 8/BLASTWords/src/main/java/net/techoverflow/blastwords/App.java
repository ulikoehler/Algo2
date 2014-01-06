package net.techoverflow.blastwords;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import org.apache.commons.io.IOUtils;
import org.python.core.PyList;
import org.python.core.PySystemState;
import org.python.util.PythonInterpreter;

/**
 * Hello world!
 *
 */
public class App {

    public static void main(String[] args) throws IOException {
        //Check correct call
        if(args.length < 5) {
            System.err.println("java -jar gruppenname blast.jar <WortLänge l> <Threshold> <Sequenz-Datei mit Sequenz in fasta-Format> <Scoring-Matrix-Datei> <Ausgabe-Datei>");
            System.exit(1);
        }
        //This is a thin wrapper around Jython.
        //See blastwords.py in the JAR for detailed sources
        PySystemState.initialize(System.getProperties(), System.getProperties(), args);
        PythonInterpreter.initialize(System.getProperties(), System.getProperties(), args);
        //Fill sys.argv
        PySystemState state = new PySystemState();
        state.argv = new PyList(Arrays.asList(args));
        System.out.println(state.argv);
        //Run the script
        PythonInterpreter interp = new PythonInterpreter(null, state);
        InputStream script = App.class.getClassLoader().getResourceAsStream("blastwords.py");
        interp.exec(IOUtils.toString(script));
    }
}
