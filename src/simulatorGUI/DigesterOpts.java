/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package simulatorGUI;

import javax.swing.JOptionPane;

/**
 *
 * @author rob
 */
public class DigesterOpts {
    private static String[] options = {"trypsin", "arg_c", "asp_n", "asp_n_ambic", "chymotrypsin", "cnbr", "lys_c", "lys_c_p", "pepsin_a", "tryp_cnbr", "tryp_chymo", "trypsin_p", "v8_de", "v8_e", "v8_e_trypsin", "v8_de_trypsin"};
    
    public static String[] getOptions(){
        return options;
    }
    
    public static Digester getDigester(String name){
        Digester result = new Digester(name,"KR","P",true);
        switch (name) {
            case "arg_c": 
                result = new Digester(name,"R","P",true);
                break;
            case "asp_n": 
                result = new Digester(name,"BD","",false);
                break;
            case "asp_n_ambic": 
                result = new Digester(name,"DE","",false);
                break;
            case "chymotrypsin":
                result = new Digester(name,"FLWY","P",true);
                break;
            case "cnbr":
                result = new Digester(name,"M","",true);
                break;
            case "lys_c":
                result = new Digester(name,"K","P",true);
                break;
            case "lys_c_p":
                result = new Digester(name,"K","",true);
                break;
            case "pepsin_a":
                result = new Digester(name,"FL","",true);
                break;
            case "tryp_cnbr":
                result = new Digester(name,"KMR","P",true);
                break;
            case "tryp_chymo":
                result = new Digester(name,"FKLRWY","P",true);
                break;
            case "trypsin_p":
                result = new Digester(name,"KR","",true);
                break;
            case "v8_de":
                result = new Digester(name,"BDEZ","P",true);
                break;
            case "v8_e":
                result = new Digester(name,"EZ","P",true);
                break;
            case "trypsin":
                result = new Digester(name,"KR","P",true);
				break;
            case "v8_e_trypsin":
                result = new Digester(name,"EKRZ","P",true);
                break;
            case "v8_de_trypsin":
                result = new Digester(name,"BDEKRZ","P",true);
                break;
            default:
                JOptionPane.showMessageDialog(null, "Unsupport digester option: " + name + ". Using trypsin.", "Error", JOptionPane.ERROR_MESSAGE);
        }
        return result;
    }
}
