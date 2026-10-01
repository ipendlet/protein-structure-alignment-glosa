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


public class AssignChemicalFeatures{
	public static void main( String[] args ) throws IOException{
		if( ( args.length == 0 ) || args[0].equalsIgnoreCase( "--help" ) ){
			System.out.println( "AssignChemicalFeatures" );
			System.out.println( "       -opts  [structure file: pdb]" );

			return;
		}
		
		FileInputStream struct_stream;
		BufferedReader struct_reader;	
		String struct_line;
		
		DataOutputStream out;

	

		// *** input file names ***
		String mStruct = args[0];

		
		String [] array;
		Atom atom;

		

		try{
			Vector<Atom> vResidue = new Vector<Atom>();
			Vector<Atom> vResidue_CF = new Vector<Atom>();
			
			array = mStruct.split( ".pdb" );
			
			out = new DataOutputStream( new FileOutputStream( new File( array[0] + "-cf.pdb" ) ) );	
			
			struct_stream = new FileInputStream( new File( mStruct ) );
			struct_reader = new BufferedReader( new InputStreamReader( struct_stream ) );	
			
			int prev_resSeq = -1000;
			int resSeq = 0;
			
			String prev_resName = "";
			String resName = "";
			
			while( ( struct_line = struct_reader.readLine()) != null ){
				if( struct_line.startsWith( "TER" ) ){
					writeChemicalFeature( prev_resName, vResidue, vResidue_CF, out );
					
					out.writeBytes( "TER" );
					
					break;
				}
				
				if( struct_line.startsWith( "ATOM" ) && !struct_line.substring( 12, 16 ).trim().startsWith( "H" ) ){
					resSeq = Integer.parseInt( struct_line.substring( 22, 26 ).trim() );
					resName = struct_line.substring( 17, 20 );
				
					if( ( resSeq != prev_resSeq ) && ( vResidue.size() != 0 ) ){
						writeChemicalFeature( prev_resName, vResidue, vResidue_CF, out );
					
						vResidue.removeAllElements();
						vResidue_CF.removeAllElements();
					}
				
					atom = new Atom();
					atom.atom_name = struct_line.substring( 12, 16 );
					atom.residue_name = struct_line.substring( 17, 20 );
					atom.chain = struct_line.substring( 21, 22 );
					atom.sequence_number = struct_line.substring( 22, 26 );
					atom.R[0] = Double.parseDouble( struct_line.substring( 30, 38 ) );
					atom.R[1] = Double.parseDouble( struct_line.substring( 38, 46 ) );
					atom.R[2] = Double.parseDouble( struct_line.substring( 46, 54 ) );		
				
					vResidue.add( atom );		
				
					prev_resSeq = resSeq;
					prev_resName = resName;
				}
			}			
			
			struct_reader.close();
			struct_stream.close();	
			
			out.close();		
		}
		catch( Exception exception ){
			exception.printStackTrace();
		}		
	}	
	
	
	
	
	
	
	public static void writeResidueInPDBFormat( Vector<Atom> vResidue, DataOutputStream out ){
		DecimalFormat format = new DecimalFormat( "####.###" );   //%8.3f
		format.setMinimumFractionDigits(3); 
		format.setMaximumFractionDigits(3); 
		
		try{ 
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				out.writeBytes( "ATOM" );
				out.writeBytes( "        " );
				out.writeBytes( vResidue.elementAt( i ).atom_name );
				out.writeBytes( " " );
				out.writeBytes( vResidue.elementAt( i ).residue_name );
				out.writeBytes( " " );
				out.writeBytes( vResidue.elementAt( i ).chain );
				out.writeBytes( vResidue.elementAt( i ).sequence_number );
				out.writeBytes( "    " );		
				
				for( int j = 0 ; j < 3 ; j++ ){
					if( vResidue.elementAt( i ).R[j] >= 0 ){
						if( vResidue.elementAt( i ).R[j] / 10 < 1 )
							out.writeBytes( "   " + format.format( vResidue.elementAt( i ).R[j] ) );
						else if( vResidue.elementAt( i ).R[j] / 10 >= 1 && vResidue.elementAt( i ).R[j] / 10 < 10 )
							out.writeBytes( "  " + format.format( vResidue.elementAt( i ).R[j] ) );
						else
							out.writeBytes( " " + format.format( vResidue.elementAt( i ).R[j] ) );
					}
					else{
						if( vResidue.elementAt( i ).R[j] / 10 > -1 )
							out.writeBytes( "  " + format.format( vResidue.elementAt( i ).R[j] ) );
						else if( vResidue.elementAt( i ).R[j] / 10 <= -1 && vResidue.elementAt( i ).R[j] / 10 > -10 )
							out.writeBytes( " " + format.format( vResidue.elementAt( i ).R[j] ) );
						else
							out.writeBytes( format.format( vResidue.elementAt( i ).R[j] ) );
					}
				}
				
				out.writeBytes( "\n" );
			}
		}
		catch( Exception exception ){
			exception.printStackTrace();
		}
	}
	
	
	
	
	
	public static void writeChemicalFeature( String residue_name, Vector<Atom> vResidue, Vector<Atom> vResidue_CF, DataOutputStream out ){
		String name = "";
		Atom atom;
		vResidue_CF.removeAllElements();
		
		if( residue_name.equalsIgnoreCase( "ARG" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "NE" ) || name.equalsIgnoreCase( "NH1" ) || name.equalsIgnoreCase( "NH2" ) ){
					atom = new Atom();
					atom.atom_name = " PC ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}
		else if( residue_name.equalsIgnoreCase( "HIS" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "ND1" ) || name.equalsIgnoreCase( "NE2" ) ){
					atom = new Atom();
					atom.atom_name = " PC ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}		
		else if( residue_name.equalsIgnoreCase( "LYS" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "NZ" ) ){
					atom = new Atom();
					atom.atom_name = " PC ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}		
		else if( residue_name.equalsIgnoreCase( "ASP" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "OD1" ) || name.equalsIgnoreCase( "OD2" ) ){
					atom = new Atom();
					atom.atom_name = " NC ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "GLU" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "OE1" ) || name.equalsIgnoreCase( "OE2" ) ){
					atom = new Atom();
					atom.atom_name = " NC ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "SER" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "OG" ) ){
					atom = new Atom();
					atom.atom_name = " OH ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}		
		else if( residue_name.equalsIgnoreCase( "THR" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "OG1" ) ){
					atom = new Atom();
					atom.atom_name = " OH ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CG2" ) ){
					atom = new Atom();
					atom.atom_name = " AL ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}
		else if( residue_name.equalsIgnoreCase( "ASN" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) || name.equalsIgnoreCase( "ND2" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) || name.equalsIgnoreCase( "OD1" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "GLN" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) || name.equalsIgnoreCase( "NE2" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) || name.equalsIgnoreCase( "OE1" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "CYS" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "SG" ) ){
					atom = new Atom();
					atom.atom_name = " AL ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "GLY" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}				
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "PRO" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "O" ) || name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}				
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "ALA" ) ){		
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CB" ) ){
					atom = new Atom();
					atom.atom_name = " AL ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}			
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "VAL" ) ){	
			double [] r = new double[3];
			r[0] = r[1] = r[2] = 0;
			
			int num_atoms_al = 0;	
			Atom al = new Atom();
					
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CB" ) || name.equalsIgnoreCase( "CG1" ) || name.equalsIgnoreCase( "CG2" ) ){
					num_atoms_al++;
					r[0] += vResidue.elementAt(i).R[0];
					r[1] += vResidue.elementAt(i).R[1];
					r[2] += vResidue.elementAt(i).R[2];
				}			
			}
			
			if( num_atoms_al == 3 ){			
				al.atom_name = " AL ";
				al.residue_name = vResidue.elementAt(0).residue_name;
				al.chain = vResidue.elementAt(0).chain;
				al.sequence_number = vResidue.elementAt(0).sequence_number;
				al.R[0] = r[0]/num_atoms_al;
				al.R[1] = r[1]/num_atoms_al;
				al.R[2] = r[2]/num_atoms_al;
				vResidue_CF.add( al );
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}		
		else if( residue_name.equalsIgnoreCase( "ILE" ) ){	
			double [] r = new double[3];
			r[0] = r[1] = r[2] = 0;
			
			int num_atoms_al = 0;	
			Atom al = new Atom();
					
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CB" ) || name.equalsIgnoreCase( "CG1" ) || name.equalsIgnoreCase( "CG2" ) || name.equalsIgnoreCase( "CD1" ) ){
					num_atoms_al++;
					r[0] += vResidue.elementAt(i).R[0];
					r[1] += vResidue.elementAt(i).R[1];
					r[2] += vResidue.elementAt(i).R[2];
				}			
			}
			
			if( num_atoms_al == 4 ){				
				al.atom_name = " AL ";
				al.residue_name = vResidue.elementAt(0).residue_name;
				al.chain = vResidue.elementAt(0).chain;
				al.sequence_number = vResidue.elementAt(0).sequence_number;
				al.R[0] = r[0]/num_atoms_al;
				al.R[1] = r[1]/num_atoms_al;
				al.R[2] = r[2]/num_atoms_al;
				vResidue_CF.add( al );
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}		
		else if( residue_name.equalsIgnoreCase( "LEU" ) ){	
			double [] r = new double[3];
			r[0] = r[1] = r[2] = 0;
			
			int num_atoms_al = 0;	
			Atom al = new Atom();
					
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CG" ) || name.equalsIgnoreCase( "CD1" ) || name.equalsIgnoreCase( "CD2" ) ){
					num_atoms_al++;
					r[0] += vResidue.elementAt(i).R[0];
					r[1] += vResidue.elementAt(i).R[1];
					r[2] += vResidue.elementAt(i).R[2];
				}			
			}
			
			if( num_atoms_al == 3 ){	
				al.atom_name = " AL ";
				al.residue_name = vResidue.elementAt(0).residue_name;
				al.chain = vResidue.elementAt(0).chain;
				al.sequence_number = vResidue.elementAt(0).sequence_number;
				al.R[0] = r[0]/num_atoms_al;
				al.R[1] = r[1]/num_atoms_al;
				al.R[2] = r[2]/num_atoms_al;
				vResidue_CF.add( al );
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "MET" ) ){						
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CE" ) ){
					atom = new Atom();
					atom.atom_name = " AL ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}		
			}			
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else if( residue_name.equalsIgnoreCase( "PHE" ) ){	
			double [] r = new double[3];
			r[0] = r[1] = r[2] = 0;
			
			int num_atoms_ar = 0;	
			Atom ar = new Atom();
					
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CG" ) || name.equalsIgnoreCase( "CD1" ) || name.equalsIgnoreCase( "CD2" ) || name.equalsIgnoreCase( "CE1" ) || name.equalsIgnoreCase( "CE2" ) || name.equalsIgnoreCase( "CZ" ) ){
					num_atoms_ar++;
					r[0] += vResidue.elementAt(i).R[0];
					r[1] += vResidue.elementAt(i).R[1];
					r[2] += vResidue.elementAt(i).R[2];
				}			
			}
			
			if( num_atoms_ar == 6 ){	
				ar.atom_name = " AR ";
				ar.residue_name = vResidue.elementAt(0).residue_name;
				ar.chain = vResidue.elementAt(0).chain;
				ar.sequence_number = vResidue.elementAt(0).sequence_number;
				ar.R[0] = r[0]/num_atoms_ar;
				ar.R[1] = r[1]/num_atoms_ar;
				ar.R[2] = r[2]/num_atoms_ar;
				vResidue_CF.add( ar );
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}		
		else if( residue_name.equalsIgnoreCase( "TYR" ) ){	
			double [] r = new double[3];
			r[0] = r[1] = r[2] = 0;
			
			int num_atoms_ar = 0;	
			Atom ar = new Atom();
					
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "OH" ) ){
					atom = new Atom();
					atom.atom_name = " OH ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}
				else if( name.equalsIgnoreCase( "CG" ) || name.equalsIgnoreCase( "CD1" ) || name.equalsIgnoreCase( "CD2" ) || name.equalsIgnoreCase( "CE1" ) || name.equalsIgnoreCase( "CE2" ) || name.equalsIgnoreCase( "CZ" ) ){
					num_atoms_ar++;
					r[0] += vResidue.elementAt(i).R[0];
					r[1] += vResidue.elementAt(i).R[1];
					r[2] += vResidue.elementAt(i).R[2];
				}			
			}
			
			if( num_atoms_ar == 6 ){	
				ar.atom_name = " AR ";
				ar.residue_name = vResidue.elementAt(0).residue_name;
				ar.chain = vResidue.elementAt(0).chain;
				ar.sequence_number = vResidue.elementAt(0).sequence_number;
				ar.R[0] = r[0]/num_atoms_ar;
				ar.R[1] = r[1]/num_atoms_ar;
				ar.R[2] = r[2]/num_atoms_ar;
				vResidue_CF.add( ar );
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}			
		else if( residue_name.equalsIgnoreCase( "TRP" ) ){	
			double [] r_1 = new double[3];			
			r_1[0] = r_1[1] = r_1[2] = 0;
			
			double [] r_2 = new double[3];			
			r_2[0] = r_2[1] = r_2[2] = 0;
			
			int num_atoms_ar_1 = 0;	
			Atom ar_1 = new Atom();
			
			int num_atoms_ar_2 = 0;	
			Atom ar_2 = new Atom();
					
			for( int i = 0 ; i < vResidue.size() ; i++ ){
				name = vResidue.elementAt(i).atom_name.trim();
			
				if( name.equalsIgnoreCase( "N" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );					
				}
				else if( name.equalsIgnoreCase( "NE1" ) ){
					atom = new Atom();
					atom.atom_name = " HD ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
					
					num_atoms_ar_2++;
					r_2[0] += vResidue.elementAt(i).R[0];
					r_2[1] += vResidue.elementAt(i).R[1];
					r_2[2] += vResidue.elementAt(i).R[2];
				}	
				else if( name.equalsIgnoreCase( "O" ) ){
					atom = new Atom();
					atom.atom_name = " HA ";
					atom.residue_name = vResidue.elementAt(i).residue_name;
					atom.chain = vResidue.elementAt(i).chain;
					atom.sequence_number = vResidue.elementAt(i).sequence_number;
					atom.R[0] = vResidue.elementAt(i).R[0];
					atom.R[1] = vResidue.elementAt(i).R[1];
					atom.R[2] = vResidue.elementAt(i).R[2];
					vResidue_CF.add( atom );
				}				
				else if( name.equalsIgnoreCase( "CE3" ) || name.equalsIgnoreCase( "CZ2" ) || name.equalsIgnoreCase( "CZ3" ) || name.equalsIgnoreCase( "CH2" ) ){
					num_atoms_ar_1++;
					r_1[0] += vResidue.elementAt(i).R[0];
					r_1[1] += vResidue.elementAt(i).R[1];
					r_1[2] += vResidue.elementAt(i).R[2];
				}	
				else if( name.equalsIgnoreCase( "CG" ) || name.equalsIgnoreCase( "CD1" ) ){
					num_atoms_ar_2++;
					r_2[0] += vResidue.elementAt(i).R[0];
					r_2[1] += vResidue.elementAt(i).R[1];
					r_2[2] += vResidue.elementAt(i).R[2];
				}
				else if( name.equalsIgnoreCase( "CD2" ) || name.equalsIgnoreCase( "CE2" ) ){
					num_atoms_ar_1++;
					r_1[0] += vResidue.elementAt(i).R[0];
					r_1[1] += vResidue.elementAt(i).R[1];
					r_1[2] += vResidue.elementAt(i).R[2];
					
					num_atoms_ar_2++;
					r_2[0] += vResidue.elementAt(i).R[0];
					r_2[1] += vResidue.elementAt(i).R[1];
					r_2[2] += vResidue.elementAt(i).R[2];
				}		
			}
			
			if( num_atoms_ar_1 == 6 ){
				ar_1.atom_name = " AR ";
				ar_1.residue_name = vResidue.elementAt(0).residue_name;
				ar_1.chain = vResidue.elementAt(0).chain;
				ar_1.sequence_number = vResidue.elementAt(0).sequence_number;
				ar_1.R[0] = r_1[0]/num_atoms_ar_1;
				ar_1.R[1] = r_1[1]/num_atoms_ar_1;
				ar_1.R[2] = r_1[2]/num_atoms_ar_1;
				vResidue_CF.add( ar_1 );
			}
			
			if( num_atoms_ar_2 == 5 ){
				ar_2.atom_name = " AR ";
				ar_2.residue_name = vResidue.elementAt(0).residue_name;
				ar_2.chain = vResidue.elementAt(0).chain;
				ar_2.sequence_number = vResidue.elementAt(0).sequence_number;
				ar_2.R[0] = r_2[0]/num_atoms_ar_2;
				ar_2.R[1] = r_2[1]/num_atoms_ar_2;
				ar_2.R[2] = r_2[2]/num_atoms_ar_2;
				vResidue_CF.add( ar_2 );
			}
			
			writeResidueInPDBFormat( vResidue_CF, out );
		}	
		else{
		}							
	}
}
