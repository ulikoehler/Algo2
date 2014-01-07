package net.techoverflow.blastwords;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
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
        //Fill sys.argv. Note that sys.argv[0] is the program name, but args[0] is not!
        PySystemState state = new PySystemState();
        LinkedList<String> argv = new LinkedList<String>(Arrays.asList(args));
        argv.addFirst("myprogram"); // = sys.argv[0]
        state.argv = new PyList(argv);
        //Run the script
        PythonInterpreter interp = new PythonInterpreter(null, state);
        InputStream script = App.class.getClassLoader().getResourceAsStream("blastwords.py");
        interp.exec(IOUtils.toString(script));
    }
}
