package infoasys.cli.pangenes;

import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import org.apache.commons.cli.*;

public class Cli {
    public Cli(String[] args) {

        Options options = new Options();

        options.addOption(Option.builder(
                "i"
        ).longOpt("input")
                .hasArg(true)
                .required()
                .desc("Input file (.faa) to process")
                .build());

        options.addOption(Option.builder(
                "k"
        ).longOpt("kvalue")
                .hasArg(true)
                .required()
                .desc("Length of the kmers used by the algorithm")
                .build());


        options.addOption(Option.builder(
                "c"
        ).longOpt("complexity")
                .hasArg(false)
                .desc("Compute the required number of operations without computing the effective network")
                .build());

        options.addOption(Option.builder(
                "o"
        ).longOpt("output")
                .hasArg(true)
                .required()
                .desc("Output file for the network")
                .build());

        options.addOption(Option.builder(
                "h"
        ).longOpt("help")
                .hasArg(false)
                .desc("Print this help message")
                .build());

        options.addOption(Option.builder(
                "j"
        ).longOpt("threads")
                .hasArg(true)
                .desc("Number of threads to use for the computation, defaults to # of processors")
                .build());

        CommandLineParser parser = new DefaultParser();

        HelpFormatter helpFormatter = new HelpFormatter();

        try {
            CommandLine cmd = parser.parse(options, args);

            if (cmd.hasOption("h")) {
                helpFormatter.printHelp("ant", options);
                System.exit(0);
            }

            InputFile = cmd.getOptionValue("i");
            OutputNetFile = cmd.getOptionValue("o");
            KValue = Integer.parseInt(cmd.getOptionValue("k"));
            OnlyOps = cmd.hasOption("c");

            if (cmd.hasOption("j")) {
                ThreadsNum = Integer.parseInt(cmd.getOptionValue("j"));
            }
            else {
                ThreadsNum = Runtime.getRuntime().availableProcessors();
            }
        }
        catch (Exception e) {
            System.out.println("Error while parsing cli arguments!");
            helpFormatter.printHelp("ant", options);
            System.exit(1);
        }
    }

    public Boolean OnlyOps;
    public String InputFile;
    public Integer KValue;
    public String OutputNetFile;
    public Integer ThreadsNum;


}
