import java.io.*;
import java.util.*;

class Atom {
    int id=0;
    String atomName;
    Atom(String an,int bt) {atomName=an; id=bt;}
}

class Bond {
    int iOtherBead=0,bondType=0,stretchingType=0;
    String otherBead;
    Bond(int ob,int bt) {
        iOtherBead=ob; bondType=bt;
    }
    Bond(String ob,int bt,int st) {
        otherBead=ob; bondType=bt; stretchingType=st;
    }
}

class Bending {
    int bendingType=0,iSecondBead=0,iThirdBead=0;
    String secondBead,thirdBead;
    Bending(int sb,int tb,int bt) {
        iSecondBead=sb; iThirdBead=tb; bendingType=bt;
    }
    Bending(String sb,String tb,int bt) {
        secondBead=sb; thirdBead=tb; bendingType=bt;
    }
}

class Torsion {
    int torsionType=0,iSecondBead=0,iThirdBead=0,iFourthBead=0;
    String secondBead,thirdBead,fourthBead;
    Torsion(int sb,int tb,int fb,int tt) {
        iSecondBead=sb; iThirdBead=tb; iFourthBead=fb; torsionType=tt;
    }
    Torsion(String sb,String tb,String fb,int tt) {
        secondBead=sb; thirdBead=tb; fourthBead=fb; torsionType=tt;
    }
}

class Bead {
    int id=0;
    String beadName;
    Set<Bond> bond=new HashSet<Bond>();
    Set<Bending> bending=new HashSet<Bending>();
    Set<Torsion> torsion=new HashSet<Torsion>();
    Bead(String name) {beadName=name;}
}

class Forcefield {
    Set<Atom> atom=new HashSet<Atom>();
    Set<Bead> bead=new HashSet<Bead>();
    Forcefield(String filename) {
        try {
            BufferedReader fin=new BufferedReader(new FileReader(filename));
            int flag=1,idx;
            boolean found;
            String line,numArrays[];
            Set<String> beadTypeSet;
            Bead nbd;
            while ((line=fin.readLine())!=null) {
                if ((idx=line.indexOf("#")) >= 0) {
                    line=line.substring(0,idx).trim();
                }
                if ((idx=line.indexOf("//")) >= 0) {
                    line=line.substring(0,idx).trim();
                }
                if ((idx=line.indexOf("!")) >= 0) {
                    line=line.substring(0,idx).trim();
                }
                numArrays=line.split(" ");
                if (numArrays.length == 1) {
                    if (numArrays[0].equalsIgnoreCase("atomtype")) flag=1;
                    else if (numArrays[0].equalsIgnoreCase("stretching")) flag=2;
                    else if (numArrays[0].equalsIgnoreCase("bending")) flag=3;
                    else if (numArrays[0].equalsIgnoreCase("torsion")) flag=4;
                    //else flag=0;
                    continue;
                }
                switch (flag) {
                case 1:
                    if (numArrays.length != 2) {System.err.println("Wrong format!"); continue;}
                    atom.add(new Atom(numArrays[0],Integer.parseInt(numArrays[1])));
                    break;
                case 2:
                    if (numArrays.length != 4) {System.err.println("Wrong format!"); continue;}
                    beadTypeSet=splitName(numArrays[0]);
                    Set<String> bondObSet=splitName(numArrays[1]);
                    Set<Integer> bondBtSet=splitNum(numArrays[2]);
                    int vibtype=Integer.parseInt(numArrays[3]);
                    for (String bdt : beadTypeSet) {
                        found=false;
                        for (Bead bd : bead) {
                            if (bd.beadName.equalsIgnoreCase(bdt)) {
                                for (String bondOb : bondObSet)
                                    for (Integer bondBt : bondBtSet)
                                        bd.bond.add(new Bond(bondOb,bondBt,vibtype));
                                found=true;
                                break;
                            }
                        }
                        if (! found) {
                            nbd=new Bead(bdt);
                            for (String bondOb : bondObSet)
                                for (Integer bondBt : bondBtSet)
                                    nbd.bond.add(new Bond(bondOb,bondBt,vibtype));
                            bead.add(nbd);
                        }
                    }
                    break;
                case 3:
                    if (numArrays.length != 4) {System.err.println("Wrong format!"); continue;}
                    beadTypeSet=splitName(numArrays[0]);
                    Set<String> bendingSbSet=splitName(numArrays[1]),bendingTbSet=splitName(numArrays[2]);
                    int bendingtype=Integer.parseInt(numArrays[3]);
                    for (String bdt : beadTypeSet) {
                        found=false;
                        for (Bead bd : bead) {
                            if (bd.beadName.equalsIgnoreCase(bdt)) {
                                for (String bendingSb : bendingSbSet)
                                    for (String bendingTb : bendingTbSet)
                                        bd.bending.add(new Bending(bendingSb,bendingTb,bendingtype));
                                found=true;
                                break;
                            }
                        }
                        if (! found) {
                            nbd=new Bead(bdt);
                            for (String bendingSb : bendingSbSet)
                                for (String bendingTb : bendingTbSet)
                                    nbd.bending.add(new Bending(bendingSb,bendingTb,bendingtype));
                            bead.add(nbd);
                        }
                    }
                    break;
                case 4:
                    if (numArrays.length != 5) {System.err.println("Wrong format!"); continue;}
                    beadTypeSet=splitName(numArrays[0]);
                    Set<String> torsionSbSet=splitName(numArrays[1]),torsionTbSet=splitName(numArrays[2]),torsionFbSet=splitName(numArrays[3]);
                    int torsiontype=Integer.parseInt(numArrays[4]);
                    for (String bdt : beadTypeSet) {
                        found=false;
                        for (Bead bd : bead) {
                            if (bd.beadName.equalsIgnoreCase(bdt)) {
                                for (String torsionSb : torsionSbSet)
                                    for (String torsionTb : torsionTbSet)
                                        for (String torsionFb : torsionFbSet)
                                            bd.torsion.add(new Torsion(torsionSb,torsionTb,torsionFb,torsiontype));
                                found=true;
                                break;
                            }
                        }
                        if (! found) {
                            nbd=new Bead(bdt);
                            for (String torsionSb : torsionSbSet)
                                for (String torsionTb : torsionTbSet)
                                    for (String torsionFb : torsionFbSet)
                                        nbd.torsion.add(new Torsion(torsionSb,torsionTb,torsionFb,torsiontype));
                            bead.add(nbd);
                        }
                    }
                    break;
                default:
                    System.out.println("Weird things do happen");
                    break;
                }
            }
            fin.close();
            System.out.println("Finished reading forcefield file");

        } catch (IOException err) {
            System.err.println("IO Exception");
        }
    }
        
    private Set<Integer> splitNum(String numStr) {
        String[] numArrays=numStr.split(",");
        Set<Integer> collect=new HashSet<Integer>();
        for (String num : numArrays) {
            collect.add(Integer.parseInt(num));
        }
        return collect;
    }

    private Set<String> splitName(String numStr) {
        String[] numArrays=numStr.split(",");
        Set<String> collect=new HashSet<String>();
        for (String num : numArrays) {
            collect.add(num);
        }
        return collect;
    } 
}
                         
class Molecule {
    private Bead bead[];
    private Forcefield forcefield;

    Molecule(String forcefieldfile,String filename) {
        int nbead=0,nbond=0,bd1,bd2,bondtype;
        String atomname;
        forcefield=new Forcefield(forcefieldfile);
        try {
            Scanner fin=new Scanner(new File(filename));
            nbead=fin.nextInt();
            bead=new Bead[nbead];
            for (int i=0;i<nbead;i++) {
                atomname=fin.next();
                bead[i]=new Bead(atomname);                        
            }
            while (true) {
                bd1=fin.nextInt(); bd2=fin.nextInt(); bondtype=fin.nextInt();
                if (bd1<1 || bd1>nbead || bd2<1 || bd2>nbead) {
                    System.err.println("Bead number "+bd1+" or "+bd2+" are out of bound -- ignore!");
                    continue;
                } else if (bondtype <1 || bondtype>3) { 
                    System.err.println("Bad edge type "+bondtype+" -- ignore!");
                    continue;
                }
                bead[bd1-1].bond.add(new Bond(bd2-1,bondtype));
                bead[bd2-1].bond.add(new Bond(bd1-1,bondtype));
                nbond++;
            }
        } catch (NoSuchElementException err) {
            if (nbond < nbead-1) {
                System.err.println("Not enough bonds!");
            } else {
                System.out.println("Finshed: nbead="+nbead+" nbond="+nbond);
                updateType();
            }
        } catch (FileNotFoundException err) {
            System.err.println("File not found");
        }
    }

    public void print(String filename) {
        try {
            PrintWriter fou=new PrintWriter(new BufferedWriter(new FileWriter(filename)));
            for (int i=0;i<bead.length;i++) {
                fou.println("! unit ntype leaderq");
                fou.printf("%-4d %-5d %d\n",i+1,bead[i].id,i+1);

                fou.println("! stretching");
                fou.println(bead[i].bond.size());
                for (Bond bd : bead[i].bond) {
                    fou.println((bd.iOtherBead+1)+" "+bd.stretchingType);
                }

                fou.println("! bending");
                fou.println(bead[i].bending.size());
                for (Bending bd : bead[i].bending) {
                    fou.println((bd.iSecondBead+1)+" "+(bd.iThirdBead+1)+" "+bd.bendingType);
                }

                fou.println("! torsion");
                fou.println(bead[i].torsion.size());
                for (Torsion ts : bead[i].torsion) {
                    fou.println((ts.iSecondBead+1)+" "+(ts.iThirdBead+1)+" "+(ts.iFourthBead+1)+" "+ts.torsionType);
                }
            }
            fou.close();
        } catch (IOException err) {
            System.err.println("IO Exception");
        }
    }

    private void updateType() {
        int path[]=new int[3];
        for (int i=0;i<bead.length;i++) {
            updateBeadType(i);
            dfs(i,0,i,path);
        }
    }

    private void updateBeadType(int cur) {
        for (Atom at : forcefield.atom) {
            if (at.atomName.equalsIgnoreCase(bead[cur].beadName)) {
                bead[cur].id=at.id;
                return;
            }
        }
    }

    private int getTorsionType(String bd0,String bd1,String bd2,String bd3) {
        for (Bead bd : forcefield.bead)
            if (bd.beadName.equalsIgnoreCase(bd0) || bd.beadName.equalsIgnoreCase(bd3)) {
                for (Torsion ts : bd.torsion)
                    if (((ts.secondBead.equalsIgnoreCase(bd1) && ts.thirdBead.equalsIgnoreCase(bd2)) || (ts.secondBead.equalsIgnoreCase(bd2) && ts.thirdBead.equalsIgnoreCase(bd1))) && ((bd.beadName.equalsIgnoreCase(bd0) && ts.fourthBead.equalsIgnoreCase(bd3)) || (bd.beadName.equalsIgnoreCase(bd3) && ts.fourthBead.equalsIgnoreCase(bd0))))
                        return ts.torsionType;
            }
        return 0;
    }

    private int getBendingType(String bd0,String bd1,String bd2) {
        for (Bead bd : forcefield.bead) {
            if (bd.beadName.equalsIgnoreCase(bd0) || bd.beadName.equalsIgnoreCase(bd2)) {
                for (Bending bdg : bd.bending)
                    if ((bd.beadName.equalsIgnoreCase(bd0) && bdg.secondBead.equalsIgnoreCase(bd1) && bdg.thirdBead.equalsIgnoreCase(bd2)) ||
                        (bd.beadName.equalsIgnoreCase(bd2) && bdg.secondBead.equalsIgnoreCase(bd1) && bdg.thirdBead.equalsIgnoreCase(bd0))) {
                        return bdg.bendingType;
                    }
            }
        }
        return 0;
    }

    private int getStretchingType(String root,String ob,int bt) {
        for (Bead bd : forcefield.bead) {
            if (bd.beadName.equalsIgnoreCase(root) || bd.beadName.equalsIgnoreCase(ob)) {
                for (Bond vib:bd.bond) if (vib.bondType==bt && ((bd.beadName.equalsIgnoreCase(root) && vib.otherBead.equalsIgnoreCase(ob)) || (bd.beadName.equalsIgnoreCase(ob) && vib.otherBead.equalsIgnoreCase(root)))) {
                        return vib.stretchingType;
                    }
            }
        }
        return 0;
    }

    private void dfs(int cur,int level,int root,int path[]) {
        switch (level) {
        case 3: bead[root].torsion.add(new Torsion(path[0],path[1],path[2],getTorsionType(bead[root].beadName,bead[path[0]].beadName,bead[path[1]].beadName,bead[path[2]].beadName))); return;
        case 2: bead[root].bending.add(new Bending(path[0],path[1],getBendingType(bead[root].beadName,bead[path[0]].beadName,bead[path[1]].beadName))); break;
        }
        for (Bond bd : bead[cur].bond) {
            if (bd.stretchingType == 0) bd.stretchingType=getStretchingType(bead[cur].beadName,bead[bd.iOtherBead].beadName,bd.bondType);
            if (((level==1) && (bd.iOtherBead==root)) || ((level==2) && ((bd.iOtherBead==path[level-2]) || (bd.iOtherBead==root)))) continue;
            path[level] = bd.iOtherBead;
            dfs(bd.iOtherBead,level+1,root,path);
        }
    }
}

public class trappebuilder {
    public static void main(final String args[]) {
        if (args.length != 3) {
            System.out.println("Usage: java trappebuilder /path/to/Forcefield /path/to/MoleculeSpec /path/to/output");
            return;
        }
        Molecule molecule=new Molecule(args[0],args[1]);
        molecule.print(args[2]);
    }
}
