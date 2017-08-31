package org.bitgun;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Main {

    private static HashMap<String, List<String>> ortologGroups = new HashMap<>();
    private static HashMap<String, List<String>> ortologFunctions = new HashMap<>();

    private static int groupId = 0;

    private static boolean createGroup(String[] tuple){
        groupId += 1;
        String id = "gTc"+groupId;
        List<String> l1 = new ArrayList<>();
        l1.add(tuple[0]);
        l1.add(tuple[1]);
        ortologGroups.put(id,l1);
        List<String> l2 = new ArrayList<>();
        l2.add(tuple[2]);
        ortologFunctions.put(id,l2);
        return true;
    }

    private static boolean addToGroup(String[] tuple){
        if(ortologGroups.size()==0){
            createGroup(tuple);
        }

        for(Map.Entry<String,List<String>> group: ortologGroups.entrySet()){
            if(group.getValue().contains(tuple[0])){
                if(!group.getValue().contains(tuple[1])){
                    group.getValue().add(tuple[1]);
                    return true;
                }
                return false;
            }else if(group.getValue().contains(tuple[1])){
                group.getValue().add(tuple[0]);
                return true;
            }
        }
        return createGroup(tuple);
    }

    public static void main(String[] args) throws FileNotFoundException {
	// write your code here

        // read line
        // check if $1 or $2 is element of a list l1;
        /// if yes; add $1 and $2 to
        /// if no; create new list and add $1,$2

        String filename = args[0];

        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line;
            while ((line = br.readLine()) != null) {
                String[] tuple = line.split("\t");
                addToGroup(tuple);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        for(Map.Entry<String,List<String>> group: ortologGroups.entrySet()){
            System.out.println("@GROUP_SIZE " + group.getValue().size());

            String sedCommand="";

            for(String gene: group.getValue()){
                String mygene = gene.replace(".","\\.");
                sedCommand += "\'s/\\<"+mygene+"/"+group.getKey()+"\\>/g\'";
                System.out.println("@MY_SED_COMMAND " + sedCommand);
                sedCommand="";

            }
        }

    }
}
