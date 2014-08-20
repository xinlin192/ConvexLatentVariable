import java.io.*;
import java.util.*;

class State{// seq1:pos seq2:(pos,k)
	
	static List<Character> bases = new ArrayList<Character>();
	static {
		bases.add('A');
		bases.add('B');
		bases.add('C');
		bases.add('D');
	}

	public int pos1;
	public int pos2;
	public char k;

	public State(int pos1, int pos2, char k){
		this.pos1 = pos1;
		this.pos2 = pos2;
		this.k = k;
	}

	public String toString(){
		return "S["+pos1+","+pos2+","+k+"]";
	}
}

class StateTrans{
	
	public State from;
	public State to;
	
	public StateTrans(State from, State to){
		this.from = from;
		this.to = to;
	}
	
	public String toString(){
		return "("+type()+","+from+", "+to+")";
	}

	public char type(){
		if( from.pos1 == to.pos1 )
			return 'D';
		else if(from.pos2 == to.pos2)
			return 'I';
		else
			return 'M';//match or mutate
	}
}

public class GenLP{
	
	static Map<String,Integer> varIndexMap = new HashMap();
	static List<String> varNameList = new ArrayList();
	
	static double lambda_trans;
	static double mut_cost;
	static double indel_cost;
	static int T_max;
	
	static double eps = 1e-2;
	static Random ran = new Random();
	
	public static void main(String[] args){
		
		if( args.length < 5 ){
			System.err.println("Usage: java GenLP [seqs] [max_T] [lambda_trans] [mut_cost] [indel_cost]");
			System.exit(0);
		}
		
		String seqsFile = args[0];
		T_max = Integer.valueOf(args[1]);
		lambda_trans = Double.valueOf(args[2]);
		mut_cost = Double.valueOf(args[3]);
		indel_cost = Double.valueOf(args[4]);
		
		List<String> seqs = new ArrayList();
		readSeqs(seqsFile, seqs);
		
		Map<Integer,Double> c_map = new HashMap(); 
		generate_c(seqs, c_map);
		System.out.println("num_var="+c_map.size());
		
		List<Double> c = new ArrayList();
		for(int i=0;i<c_map.size();i++){
			c.add(c_map.get(i));
		}
		
		List<List<Pair<Integer,Double>>> Aeq = new ArrayList();
		List<Double> beq = new ArrayList();
		generate_Aeq_beq(seqs, Aeq, beq);

		List<List<Pair<Integer,Double>>> A = new ArrayList();
		List<Double> b = new ArrayList();
		generate_A_b(seqs, A, b);
		
		List<Double> lb = new ArrayList();
		for(int i=0;i<varNameList.size();i++){
			lb.add(0.0);
		}

		writeLP(c,A,b,Aeq,beq,lb);
		writeVarName("varName",varNameList);
	}
	
	static void writeVarName(String fname, List<String> names){
		try{
			BufferedWriter bufw = new BufferedWriter(new FileWriter(fname));
			for(String name: names){
				bufw.write(name+"\n");
			}
			bufw.close();

		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

	static void writeLP(List<Double> c, List<List<Pair<Integer,Double>>> A, List<Double> b, List<List<Pair<Integer,Double>>> Aeq, List<Double> beq, List<Double> lb){
		
		writeVect("c",c);
		writeMat("A",A);
		writeMat("Aeq",Aeq);
		writeVect("b",b);
		writeVect("beq",beq);
		writeVect("lb",lb);
	}
	
	static void writeMat(String fname, List<List<Pair<Integer,Double>>> A){
		
		try{
			BufferedWriter bufw = new BufferedWriter(new FileWriter(fname));
			bufw.write(A.size()+"\t"+varNameList.size()+"\t"+0+"\n");
			for(int i=0;i<A.size();i++){
				for(int j=0;j<A.get(i).size();j++){
					bufw.write((i+1)+"\t"+(A.get(i).get(j).first+1)+"\t"+A.get(i).get(j).second+"\n");
				}
			}
			bufw.close();

		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}

	static void writeVect(String fname, List<Double> b){
		try{
			BufferedWriter bufw = new BufferedWriter(new FileWriter(fname));
			for(Double v: b){
				bufw.write(v+"\n");
			}
			bufw.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	static void generate_A_b(List<String> seqs, List<List<Pair<Integer,Double>>> A, List<Double> b){
		
		State sf, st;
		StateTrans trans;
		int index;
		String seq;
		List<Pair<Integer,Double>> a;
		for(Integer n=0;n<seqs.size();n++){
			seq = seqs.get(n);
			for(int i=0;i<seq.length();i++){
				for(int j=0;j<T_max;j++){
					for(Character k: State.bases){
						
						
						sf = new State(i,j,k);
						//deletion
						if(j+1<T_max){
							for(Character k2: State.bases){
								st = new State(i,j+1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString()); //n, sf, st
								
								a = new ArrayList();
								a.add(new Pair(index,1.0));
								index = varIndex("xi["+j+","+k+","+k2+"]");
								a.add(new Pair(index,-1.0));

								A.add(a);
								b.add(0.0);
							}
						}
						
						//mutation or match
						if(i+1<seq.length() && j+1<T_max){
							for(Character k2: State.bases){
								st = new State(i+1,j+1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString()); //n, sf, st
							
								a = new ArrayList();
								a.add(new Pair(index,1.0));
								index = varIndex("xi["+j+","+k+","+k2+"]");
								a.add(new Pair(index,-1.0));
								
								A.add(a);
								b.add(0.0);
							}
						}
						if( i==seq.length()-1 ){//consider whether j trans to 'END'
							st = new State(i+1,j+1,'$');
							trans = new StateTrans(sf,st);
							index = varIndex(n+trans.toString());
							
							a = new ArrayList();
							a.add(new Pair(index,1.0));
							index = varIndex("xi["+j+","+k+",$]");
							a.add(new Pair(index,-1.0));
							
							A.add(a);
							b.add(0.0);
						}
					}
					
				}
			}
		}
	}
	
	static void generate_Aeq_beq(List<String> seqs, List<List<Pair<Integer,Double>>> Aeq, List<Double> beq){
		
		//flows from (1,1) sum to 1
		List<Pair<Integer,Double>> a ;
		State sf, st;
		StateTrans trans;
		int index;
		for(Integer n=0;n<seqs.size();n++){
			a = new ArrayList();
			for(Character k: State.bases){
				//insertioin
				sf = new State(0,0,k);
				st = new State(1,0,k);
				trans = new StateTrans(sf,st);
				index = varIndex(n + trans.toString());
				a.add(new Pair(index,1.0));

				//deletion
				for(Character k2: State.bases){
					st = new State(0,1,k2);
					trans = new StateTrans(sf,st);
					index = varIndex(n + trans.toString());
					a.add(new Pair(index,1.0));
				}

				//mutation, match
				for(Character k2: State.bases){
					st = new State(1,1,k2);
					trans = new StateTrans(sf,st);
					index = varIndex(n + trans.toString());
					a.add(new Pair(index,1.0));
				}
			}
			Aeq.add(a);
			beq.add(1.0);
		}
		String seq;
		for(Integer n=0;n<seqs.size();n++){
			seq = seqs.get(n);
			//go through the |seq| by T_max matrix of "states"
			for(int i=0;i<seq.length();i++){
				for(int j=0;j<T_max;j++){
					
					if(i==0 && j==0)
						continue;
					
					for(Character k: State.bases){
						
						a = new ArrayList();
						
						// to this state
						st = new State(i,j,k);
						////insertion
						if( i-1 >= 0 ){
							sf = new State(i-1,j,k);
							trans = new StateTrans(sf,st);
							index = varIndex(n+trans.toString());
							a.add(new Pair(index,1.0));
						}
						////deletion
						if( j-1 >= 0 ){
							for(Character k2: State.bases){
								sf = new State(i,j-1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString());
								a.add(new Pair(index,1.0));
							}
						}
						////mutation/match
						if( i-1 >=0 && j-1 >= 0 )
							for(Character k2: State.bases){
								sf = new State(i-1,j-1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString());
								a.add(new Pair(index,1.0));
							}

						// from this state
						sf = new State(i,j,k);
						////insertion
						if( i+1 < seq.length() ){
							st = new State(i+1,j,k);
							trans = new StateTrans(sf,st);
							index = varIndex(n+trans.toString());
							a.add(new Pair(index,-1.0));
						}
						////deletion
						if( j+1 < T_max ){
							for(Character k2: State.bases){
								st = new State(i,j+1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString());
								a.add(new Pair(index,-1.0));
							}
						}
						////mutate/match
						if( i+1 < seq.length() && j+1 < T_max ){
							for(Character k2: State.bases){
								st = new State(i+1,j+1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString());
								a.add(new Pair(index,-1.0));
							}
						}
						if( i==seq.length()-1 ){
							st = new State(i+1,j+1,'$');
							trans = new StateTrans(sf,st);
							index = varIndex(n+trans.toString());
							a.add(new Pair(index,-1.0));
						}

						Aeq.add(a);
						beq.add(0.0);
					}
				}
			}
		}
		
	}

	static void generate_c(List<String> seqs, Map<Integer,Double> c){
		
		String seq;
		State sf, st;
		StateTrans trans;
		int index;
		//c_{n,sf,st}
		for(Integer n=0;n<seqs.size();n++){
			seq = seqs.get(n);
			//go through the |seq| by T_max matrix of "states"
			for(int i=0;i<seq.length();i++){
				for(int j=0;j<T_max;j++){
					for(Character k: State.bases){
						
						sf = new State(i,j,k);
						
						//insertion
						if(i+1<seq.length()){
							st = new State(i+1,j,k);
							trans = new StateTrans(sf,st);
							index = varIndex(n+trans.toString()); //n, sf, st
							c.put(index,indel_cost+noise());
						}

						//deletion
						if(j+1<T_max){
							for(Character k2: State.bases){
								st = new State(i,j+1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString()); //n, sf, st
								c.put(index,indel_cost+noise());
							}
						}

						//mutation or match
						double cost = (seq.charAt(i)==k)?0.0:mut_cost;
						if(i+1<seq.length() && j+1<T_max){
							for(Character k2: State.bases){
								st = new State(i+1,j+1,k2);
								trans = new StateTrans(sf,st);
								index = varIndex(n+trans.toString()); //n, sf, st
								c.put(index,cost+noise());
							}
						}
						if( i==seq.length()-1 ){//consider whether j trans to 'END'
							st = new State(i+1,j+1,'$');
							trans = new StateTrans(sf,st);
							index = varIndex(n+trans.toString());
							c.put(index,cost+noise());
						}
					}
				}
			}
		}

		//c_{xi_{t,i,j}}
		for(int t=0;t<T_max;t++){
			for( Character i: State.bases ){
				if(t+1<T_max){
					for(Character j: State.bases){
						index = varIndex("xi["+t+","+i+","+j+"]");
						c.put(index,lambda_trans+noise());
					}
				}
				index = varIndex("xi["+t+","+i+",$]");
				c.put(index,lambda_trans+noise());
			}
		}
	}

	static double noise(){
		return eps * ran.nextGaussian();
	}

	static Integer varIndex(String varName){
		
		Integer ind =  varIndexMap.get(varName);
		if( ind != null )
			return ind;
		else{
			ind = varNameList.size();
			varNameList.add(varName);
			varIndexMap.put(varName,ind);
			return ind;
		}
	}

	static void readSeqs(String file, List<String> seqs){
		
		try{
			BufferedReader bufr = new BufferedReader(new FileReader(file));

			String line;
			while( (line=bufr.readLine()) != null ){
				seqs.add(line);
			}
			bufr.close();
			
		}catch(Exception e){
			e.printStackTrace();
			System.exit(0);
		}
	}
}
