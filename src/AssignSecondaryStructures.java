import java.io.*;
import java.util.Vector;
import java.text.DecimalFormat;


class Atom{
	String atom_name;        
	String residue_name;     
	String chain;
	String sequence_number; 
	char ss;

	double[] R;

	public Atom(){
		ss = 'X';
		
		R = new double[3];	

		R[0] = 0.0;
		R[1] = 0.0;
		R[2] = 0.0;
	}
}


public class AssignSecondaryStructures{
	static double lambda_A[] = { 0, 0, 5.45, 5.18, 6.37 };
	static double lambda_B[] = { 0, 0, 6.1, 10.4, 13 };
	static double delta_A = 2.1;
	static double delta_B = 1.42;
	static char [] secondaryStructure = new char[7000];
	
	public static void main( String[] args ) throws IOException{
		if( ( args.length == 0 ) || args[0].equalsIgnoreCase( "--help" ) ){
			System.out.println( "AssignSecondaryStructures" );
			System.out.println( "       -opts  [local  structure file: pdb]" );
			System.out.println( "              [global structure file: pdb]" );

			return;
		}
		
		FileInputStream struct_stream;
		BufferedReader struct_reader;	
		String struct_line;
		
		DataOutputStream out;
		
		String mLocalStruct = args[0];
		String mGlobalStruct = args[1];

		
		String [] array;
		Atom atom;

		

		try{
			// Read input structure
			Vector<Atom> vLocalStruct = new Vector<Atom>();
			
			struct_stream = new FileInputStream( new File( mLocalStruct ) );
			struct_reader = new BufferedReader( new InputStreamReader( struct_stream ) );	
			
			while( ( struct_line = struct_reader.readLine()) != null ){
				if( struct_line.startsWith( "TER" ) )
					break;
				
				if( struct_line.startsWith ( "ATOM" ) && struct_line.substring( 13, 15 ).equalsIgnoreCase( "CA" ) ){
					atom = new Atom();
					atom.atom_name = struct_line.substring( 12, 16 );
					atom.residue_name = struct_line.substring( 17, 20 );
					atom.chain = struct_line.substring( 21, 22 );
					atom.sequence_number = struct_line.substring( 22, 26 );
					atom.R[0] = Double.parseDouble( struct_line.substring( 30, 38 ) );
					atom.R[1] = Double.parseDouble( struct_line.substring( 38, 46 ) );
					atom.R[2] = Double.parseDouble( struct_line.substring( 46, 54 ) );	
					
					vLocalStruct.add( atom );		
				}
			}			
			
			struct_reader.close();
			struct_stream.close();	
			
			
			
			// Read input global structure
			Vector<Atom> vGlobalStruct = new Vector<Atom>();
			
			struct_stream = new FileInputStream( new File( mGlobalStruct ) );
			struct_reader = new BufferedReader( new InputStreamReader( struct_stream ) );	
			
			while( ( struct_line = struct_reader.readLine()) != null ){
				if( struct_line.startsWith( "TER" ) )
					break;
				
				if( struct_line.startsWith ( "ATOM" ) && struct_line.substring( 13, 15 ).equalsIgnoreCase( "CA" ) ){
					atom = new Atom();
					atom.atom_name = struct_line.substring( 12, 16 );
					atom.residue_name = struct_line.substring( 17, 20 );
					atom.chain = struct_line.substring( 21, 22 );
					atom.sequence_number = struct_line.substring( 22, 26 );
					atom.R[0] = Double.parseDouble( struct_line.substring( 30, 38 ) );
					atom.R[1] = Double.parseDouble( struct_line.substring( 38, 46 ) );
					atom.R[2] = Double.parseDouble( struct_line.substring( 46, 54 ) );	
					
					vGlobalStruct.add( atom );		
				}
			}			
			
			struct_reader.close();
			struct_stream.close();	
			
			
			
			setSecondaryStructure( vGlobalStruct, vLocalStruct );
			
			
			
			array = mLocalStruct.split( ".pdb" );			
			out = new DataOutputStream( new FileOutputStream( new File( array[0] + "-ss.pdb" ) ) );	
			
			for( int i = 0 ; i < vLocalStruct.size() ; i++ ){
				out.writeBytes( vLocalStruct.elementAt(i).sequence_number + " " + vLocalStruct.elementAt(i).residue_name + " " + " CA " + vLocalStruct.elementAt(i).ss + "\n" );
			}
			
			out.writeBytes( "TER" );
			out.close();		
		}
		catch( Exception exception ){
			exception.printStackTrace();
		}		
	}	
	
	
	
	
	
	
	public static void resetSecondaryStructureArray(){
		for( int i = 0 ; i < 3000; i++ )
			secondaryStructure[i] = 'X';
	}
	
	
	
	
	
	public static void setSecondaryStructure( Vector< Atom > vGlobalStruct, Vector< Atom > vLocalStruct ){
		resetSecondaryStructureArray();
	
		int n = vGlobalStruct.size() - 4;
		double dist;
		double delta;
	
		int isAlpha = 1;
		int isBeta = 1;
	
		for( int i = 2 ; i < n ; i++ ){
			isAlpha = 1;
			isBeta = 1;
	
			for( int j = i - 2 ; j <= i ; j++ ){
				for( int k = 2 ; k <= 4 ; k++ ){
					dist = Math.sqrt( Math.pow( vGlobalStruct.elementAt(j).R[0] - vGlobalStruct.elementAt(j+k).R[0], 2 ) + Math.pow( vGlobalStruct.elementAt(j).R[1] - vGlobalStruct.elementAt(j+k).R[1], 2 ) + Math.pow( vGlobalStruct.elementAt(j).R[2] - vGlobalStruct.elementAt(j+k).R[2], 2 ) );
				
					if( Math.abs( dist - lambda_A[k] ) >= delta_A )
						isAlpha = 0;
					
					if( Math.abs( dist - lambda_B[k] ) >= delta_B )
						isBeta = 0;				
				}
			}
		
			if( Integer.parseInt( vGlobalStruct.elementAt(i).sequence_number.trim() ) >= 0 ){		
				if( isAlpha == 1 )
					secondaryStructure[Integer.parseInt( vGlobalStruct.elementAt(i).sequence_number.trim() )] = 'A';
				else if( isBeta == 1 )
					secondaryStructure[Integer.parseInt( vGlobalStruct.elementAt(i).sequence_number.trim() )] = 'B';
				else 
					secondaryStructure[Integer.parseInt( vGlobalStruct.elementAt(i).sequence_number.trim() )] = 'C';
			}
		}
	
	
		n = vLocalStruct.size();
	
		for( int i = 0 ; i < n ; i++ ){
			if( Integer.parseInt( vLocalStruct.elementAt(i).sequence_number.trim() ) >= 0 )			
				vLocalStruct.elementAt(i).ss = secondaryStructure[Integer.parseInt( vLocalStruct.elementAt(i).sequence_number.trim() )];
			else
				vLocalStruct.elementAt(i).ss = 'C';
		}
	}
}
