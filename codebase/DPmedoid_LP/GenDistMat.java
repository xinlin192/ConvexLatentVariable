import java.io.*;
import java.util.*;

public class GenDistMat{
	
	public static void main(String[] args){
		
		if( args.length < 1 ){
			System.err.println("Usage: java GenDistMat [data]");
			System.exit(0);
		}

		String datafile = args[0];
		List<Instance> data = readData(datafile);
		System.err.println("data size=" + data.size());
		List<double[]> distMat = compute_dist_mat(data);
		
		normalize(distMat);

		writeToFile(distMat, datafile+".dist");
	}
	
	static void normalize(List<double[]> distMat){
		
		double maxVal = maxVal(distMat);
		for(int i=0;i<distMat.size();i++){
			double[] dist_arr = distMat.get(i);
			for(int j=0;j<dist_arr.length;j++){
				dist_arr[j] /= maxVal;
			}
		}
	}

	static double maxVal(List<double[]> distMat){
		
		double maxVal = -1e300;
		for(double[] dist_arr: distMat){
			for(int j=0;j<dist_arr.length;j++){
				if( dist_arr[j] > maxVal) 
					maxVal = dist_arr[j];
			}
		}
		return maxVal;
	}

	static void writeToFile(List<double[]> distMat, String fname){
		
		try{
			BufferedWriter bufw = new BufferedWriter(new FileWriter(fname));
			for(int i=0;i<distMat.size();i++){
				for(int j=0;j<distMat.get(i).length;j++)
					bufw.write(distMat.get(i)[j]+" ");
				bufw.newLine();
			}
			bufw.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

	static List<double[]> compute_dist_mat(List<Instance> data){
		
		int N = data.size();
		List<double[]> mat = new ArrayList();
		for(int i=0;i<N;i++){
			double[] dist_arr = new double[N];
			for(int j=0;j<N;j++){
				dist_arr[j] = dist(data.get(i),data.get(j));
			}
			mat.add(dist_arr);
		}
		return mat;
	}

	static double dist(Instance a, Instance b){
		
		//Euclidian distance
		int i=0,j=0;
		double sum = 0.0;
		Pair<Integer,Double> kv_a;
		Pair<Integer,Double> kv_b;
		while( i<a.size() && j<b.size() ){
			kv_a = a.get(i);
			kv_b = b.get(j);
			if( kv_a.key.equals(kv_b.key) ){
				double dif = (kv_a.value - kv_b.value);
				sum += dif*dif;
				i++; 
				j++;
			}
			else if( kv_a.key < kv_b.key ){
				double dif = kv_a.value;
				sum += dif*dif;
				i++;
			}
			else if( kv_a.key > kv_b.key ){
				double dif = -kv_b.value;
				sum += dif*dif;
				j++;
			}
		}

		if( i < a.size() ){
			double dif;
			for(;i<a.size();i++){
				dif = a.get(i).value;
				sum += dif*dif;
			}
		}else if( j < b.size() ){
			double dif;
			for(;j<b.size();j++){
				dif = -b.get(i).value;
				sum += dif*dif;
			}
		}
		
		return Math.sqrt(sum);
		//return sum;
	}

	static List<Instance> readData(String fname){
		
		List<Instance> data = new ArrayList();
		try{
			BufferedReader bufr = new BufferedReader(new FileReader(fname));
			
			String line;
			String[] tokens;
			String[] kv;
			while( (line=bufr.readLine()) != null ){
				tokens = line.trim().split(" ");
				Instance ins = new Instance();
				for(int i=1;i<tokens.length;i++){
					if( tokens[i].length() < 1 )
						continue;
					kv = tokens[i].split(":");
					ins.add( new Pair(Integer.valueOf(kv[0]), Double.valueOf(kv[1])) );
				}
				data.add(ins);
			}

		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
			
		return data;
	}
}


class Instance extends ArrayList<Pair<Integer,Double>>{}
