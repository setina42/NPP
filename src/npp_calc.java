import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.io.*;
import java.util.concurrent.TimeUnit;

import org.w3c.dom.Document;
import org.w3c.dom.*;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;


import org.apache.commons.text.similarity.LevenshteinDistance;

public class npp_calc {

	//Initial Params
	//public static double BlastSimThreshhold = 70;
	//public static int NPPLengthDiffThreshhold = 30;
	//public static int NPPLengthThreshhold = 500;
	//public static String GeneInputFile = "/data/Copci_AmutBmut1_GeneModels_FrozenGeneCatalog_20160912_aa.fasta";
	//public static String OutputFile = "output_Sequences"; 	
	
	
	// Params for short NPPs
	public static double BlastSimThreshhold = 40;
	public static int NPPLengthDiffThreshhold = 8;
	public static int NPPLengthThreshhold = 20;
	public static String GeneInputFile = "Copci_AmutBmut1_GeneModels_FrozenGeneCatalog_20160912_aa.fasta";
	public static String OutputFile = "output_ShortSequences"; 

    public static void main(String[] args) throws IOException, InterruptedException {
        new npp_calc().go();
    }

    private void go() throws IOException, InterruptedException {

        /** Read the input data */
		BufferedReader in = new BufferedReader(new FileReader(GeneInputFile));

        String str;
        StringBuilder sb = new StringBuilder();
        while ((str = in.readLine()) != null)
            sb.append(str);
        in.close();

        String data = sb.toString();


        //split data into genes (remove all the other stuff
        String[] genes = data.split(">jgi");
        FileWriter writer = new FileWriter(OutputFile + "" + ".csv");
        writer.append("Blast Similars;Equals;same length;smart substring;substring;levenstein;used ld-distance;sequence;gene name;completSequence\n");
        int progress = 0;
        for (String x : genes) {

            progress++;
            System.out.println(progress + " of " + genes.length + " genes");

            if (x == null || x.isEmpty()) {
                continue;
            }
            String geneName = x.substring(0, 25);
            if (x.length() > 800 ){
                continue;
            }
            //writer.append("new gene\n");

            /** Split the data */
            final int IGNORE_IF_SMALLER_THAN = 3;
            String[] cleavingSites = x.split("(KK|KR|RR|RK)");
            cleavingSites[0] = cleavingSites[0].substring(17).replaceAll("[a0-9]", "");
            Set<String> set = new LinkedHashSet<>();

            for (String s : cleavingSites) {
                if (s.length() >= IGNORE_IF_SMALLER_THAN)
                    set.add(s);
            }
            //no KK, KR... found
            if (set.size() == 0) {
                continue;
            }

            // System.out.println("Found " + set.size() + " different Strings");

            List<MyString> uniqueList = new ArrayList<>();

            /** Get rid of duplicates using Set*/
            for (String s : set)
                uniqueList.add(new MyString(s));

            /** compare each String to all from this gene using the equality functions */
            int count = 0;
            long lastPercent = 0;


            long startMillis = System.currentTimeMillis();
            for (MyString unique : uniqueList) {
                count++;
                if (unique.s.length() > 20) {
                    continue;
                }
                //int BlastSimilars = 0;
                long percent = 100 * count / uniqueList.size();
//                if (Math.floor(percent) != Math.floor(lastPercent)) {
//                    long currentMillis = System.currentTimeMillis();
//                    System.out.println("Processed " + percent + "% (Remaining time: " + (currentMillis - startMillis) / 1000.0 / 60.0 / percent * (100.0 - percent) + "min)");
//                    lastPercent = percent;
//                }
                //System.out.println("### Processing unique: " + unique.s);
                for (String cleavingSite : cleavingSites) {
                    //System.out.println("Processing cleaving: " +cleavingSite);
                    unique.evalSimilarity1(cleavingSite);
                    unique.evalSimilarity2(cleavingSite);
                    //unique.evalSimilarity3(cleavingSite);
                    //unique.evalSimilarity4(cleavingSite);
                    //unique.evalSimilarity5(cleavingSite);

                    //make some constraints about when a blast search makes sense, similar length:
                    // neuropeptide: assume similar length between individual string

                    if ((unique.s.length() - cleavingSite.length()) < NPPLengthDiffThreshhold) {
                        unique.evalBlast(cleavingSite);
                    }

                    //blast search
                    //output both string to compare to a txt.file for further use in Blast+ standalone toolkit
//
//                        FileWriter site1 = new FileWriter("site1" + "" + ".txt");
//                        site1.append(cleavingSite);
//                        site1.flush();
//                        site1.close();
//                        FileWriter site2 = new FileWriter("site2" + "" + ".txt");
//                        site2.append(unique.s);
//                        site2.flush();
//                        site2.close();
//					Runtime rt = Runtime.getRuntime();
//					String[] commands = {"blastp","-query", "site1.txt", "-subject", "site2.txt",  "-outfmt", "6 pident"};
//                    //String[] commands = {"blastp","-query", "site1.txt", "-subject", "site2.txt", "-out", "alignment.xml", "-outfmt", "5"};
//					//String[] commands = {"blastp","-help"};
//					Process proc = rt.exec(commands);
//
//					BufferedReader stdInput = new BufferedReader(new
//							InputStreamReader(proc.getInputStream()));
//
//					BufferedReader stdError = new BufferedReader(new
//							InputStreamReader(proc.getErrorStream()));
//
//					boolean HitFound = false;
//   				read the output from the command
//					System.out.println("Here is the standard output of the command:\n");
//					String s = null;
//					while ((s = stdInput.readLine()) != null) {
//						System.out.println(s);
//						//if(s.split("No hits found").length > 1){
//						    HitFound = true;
//						    int a = 3;
//						//}
//					}
//
// 				read any errors from the attempted command
//					System.out.println("Here is the standard error of the command (if any):\n");
//					while ((s = stdError.readLine()) != null) {
//						//System.out.println(s);
//					}
//					if(HitFound=true){BlastSimilars++;}
////

                }

                if (unique.blastsimilars > 1 || unique.equal > 1) {
                    try {
                        // Create File
                        writer.append(unique.blastsimilars + ";" + unique.equal + ";" + unique.length + ";" + unique.subs + ";" + unique.sub + ";" + unique.levenshtein + ";;" + unique.s + ";" + geneName + ";"+ x+ "\n");
                    } catch (Exception e) {
                        // PrintError
                        System.out.println(e.getMessage());
                    }
                }
                //writer.append(unique.equal + ";" + unique.length + ";" + unique.subs + ";" + unique.sub + ";" + unique.levenshtein + ";" + unique.ld.getThreshold() + ";" + unique.s + "\n");

            }

        }
        writer.flush();
        writer.close();
        System.out.println("Done :))");
    }
    
    class MyString {
        String s;
        Integer equal = 0;
        Integer length = 0;
        Integer subs = 0;
        Integer sub = 0;
        Integer levenshtein = 0;
        Integer blastsimilars = 0;

        LevenshteinDistance ld;

        public MyString(String s) {
            this.s = s;
        }

        void evalSimilarity1(String other) {
            if (s.equals(other))
                equal++;
        }

        void evalSimilarity2(String other) {
            if (Math.abs(s.length() - other.length()) < 5)
                length++;
        }

        void evalSimilarity3(String other) {
            if (other.equals(s))
                return;

            if (other.equals(""))
                return;

            int io = 0;
            for (int i = 0; i < s.length(); i++) {
                if (io > other.length() - 1)
                    return;
                while (!(s.charAt(i) == other.charAt(io))) {
                    io++;
                    if (io > other.length() - 1)
                        return;
                }
                io++;
            }
            subs++;

        }

        void evalSimilarity4(String other) {
            if (other.equals(s))
                return;
            if (other.contains(s))
                sub++;

        }

        void evalBlast(String other) throws IOException {


            try {
                FileWriter site1 = new FileWriter("site1" + "" + ".txt");
                site1.append(other);
                site1.flush();
                site1.close();
            } catch (FileNotFoundException e) {

                try {
                    Thread.sleep(2000);
                    FileWriter site1 = new FileWriter("site1" + "" + ".txt");
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();
                }
                FileWriter site1 = new FileWriter("site1" + "" + ".txt");
                site1.append(other);
                site1.flush();
                site1.close();
            }


            try {
                FileWriter site2 = new FileWriter("site2" + "" + ".txt");
                site2.append(s);
                site2.flush();
                site2.close();
            } catch (FileNotFoundException e) {

                try {
                    Thread.sleep(2000);
                    FileWriter site2 = new FileWriter("site2" + "" + ".txt");
                } catch (InterruptedException ex) {
                    Thread.currentThread().interrupt();

                    FileWriter site1 = new FileWriter("site1" + "" + ".txt");
                    site1.append(other);
                    site1.flush();
                    site1.close();
                }
            }


            Runtime rt = Runtime.getRuntime();
            String[] commands = {"blastp", "-query", "site1.txt", "-subject", "site2.txt", "-outfmt", "6 pident"};
            //String[] commands = {"blastp","-query", "site1.txt", "-subject", "site2.txt", "-out", "alignment.xml", "-outfmt", "5"};
            //String[] commands = {"blastp","-help"};
            Process proc = rt.exec(commands);

            BufferedReader stdInput = new BufferedReader(new
                    InputStreamReader(proc.getInputStream()));

            BufferedReader stdError = new BufferedReader(new
                    InputStreamReader(proc.getErrorStream()));

            boolean HitFound = false;
            // read the output from the command
            //System.out.println("Here is the standard output of the command:\n");
            String s = null;
            while ((s = stdInput.readLine()) != null) {
                //System.out.println(s);
                //if(s.split("No hits found").length > 1){
                if ((Double.parseDouble(s)) > BlastSimThreshhold) {
                    blastsimilars++;
                }

                //}
            }

            // read any errors from the attempted command
            //System.out.println("Here is the standard error of the command (if any):\n");
            while ((s = stdError.readLine()) != null) {
                //System.out.println(s);
            }
        }


        void evalSimilarity5(String other) {
            if (other.equals(s))
                return;

            if (ld == null) {
                int highestDistance = (int) Math.max(2, (0.05 * s.length()));
                ld = new LevenshteinDistance(highestDistance);
            }

            if (!ld.apply(s, other).equals(-1))
                levenshtein++;
        }
    }

}