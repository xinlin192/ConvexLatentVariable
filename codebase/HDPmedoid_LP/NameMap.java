import java.io.*;
import java.util.*;

class NameMap{
	public static void main(String[] args){
		if( args.length < 2 ){
			System.err.println("java NameMap [sol] [varName]");
			System.exit(0);
		}

		String sol_file = args[0];
		String name_file = args[1];
		
		List<String> nameList = new ArrayList();
		readList(name_file, nameList);
		
		try{
			BufferedReader bufr =new BufferedReader(new FileReader(sol_file));
			BufferedWriter bufw =new BufferedWriter(new FileWriter(sol_file+".name"));
			String line;
			String[] tokens;
			Integer index;
			while( (line=bufr.readLine()) != null ){
				tokens = line.split("\t");
				if(tokens.length==2){
					index = Integer.valueOf(tokens[0])-1;
					bufw.write(nameList.get(index)+"\t"+tokens[1]+"\n");
				}else if(tokens.length==3){
					index = Integer.valueOf(tokens[1]);
					bufw.write(tokens[0]+"\t"+nameList.get(index)+"\t"+tokens[2]+"\n");
				}
			}
			bufw.close();
			bufr.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

	static void readList(String fname, List<String> list){
		
		try{
			BufferedReader bufr = new BufferedReader(new FileReader(fname));
			String line;
			while( (line=bufr.readLine()) != null ){
				list.add(line);
			}
			bufr.close();

		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}
}
