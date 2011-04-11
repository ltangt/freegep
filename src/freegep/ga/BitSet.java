package freegep.ga;

import java.util.Arrays;

public class BitSet {

	boolean[] bits;
	
	public BitSet(int length){
		bits = new boolean[length];
		Arrays.fill(bits, false);
	}
	
	public BitSet clone(){
		boolean[] newBits = Arrays.copyOf(bits, bits.length);
		BitSet newSet = new BitSet(bits.length);
		newSet.bits = newBits;
		return newSet;
	}
	
	public void set(int idx, boolean value){
		bits[idx] = value;
	}
	
	public void set(int idx){
		set(idx, true);
	}
	
	public boolean get(int idx){
		return bits[idx];
	}
	
	public int length(){
		return bits.length;
	}
	
	public void or(BitSet bitSet){
		boolean[] otherBits = bitSet.bits;
		for(int i = 0; i < bits.length; ++i){
			if(true == otherBits[i]){
				bits[i] = true;
			}
		}
	}
	
	public void and(BitSet bitSet){
		boolean[] otherBits = bitSet.bits;
		for(int i = 0; i < bits.length; ++i){
			if(true == otherBits[i] && true == bits[i]){
				bits[i] = true;
			}
			else{
				bits[i] = false;
			}
		}
	}
	
	public boolean equals(Object o){
		if(o instanceof BitSet){
			BitSet oBitSet = (BitSet)o;
			if(oBitSet.length() == bits.length){
				boolean[] oBits = oBitSet.bits;
				for(int i = 0; i < oBits.length; ++i){
					if(oBits[i] != bits[i]){
						return false;
					}
				}
				return true;
			}
		}
		
		return false;
	}
	
	public String toString(){
		StringBuffer output = new StringBuffer();
		for(int i = 0; i < bits.length; ++i){
			if(true == bits[i]){
				output.append('1');
			}
			else{
				output.append('0');
			}
		}
		return output.toString();
	}
	
	public static void main(String[] args){
		int len = 10;
		BitSet b1 = new BitSet(len);
		BitSet b2 = new BitSet(len);
		BitSet b3 = new BitSet(len);
		b3.set(3);
		
		System.out.println(b1.equals(b2));
		System.out.println(b2.equals(b3));
		System.out.println(b2);
		b2.set(3);
		System.out.println(b1.equals(b2));
		System.out.println(b2.equals(b3));
		System.out.println(b2);
	}
	
}
