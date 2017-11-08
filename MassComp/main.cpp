/*
MassComp: mass spectrometry data compression

Ruochen Yang
email: rcyang624@126.com
*/

#include "head.h"
#ifndef INT_RESEARCH_RANGE
#define INT_SEARCH_RANGE 15
#endif

#ifndef MAX_PATH
#define MAX_PATH 200
#endif

// interger implementation of arithmetic coding: 
// M is the lower,upper,tag's binary length
#ifndef M
//#define M 16
#define M 32
#endif

// using macro for checking error results
#ifndef XMLCheckResult
#define XMLCheckResult(a_eResult) if (a_eResult != XML_SUCCESS) { printf("Error: %i\n", a_eResult); return a_eResult; }
#endif




void Comp::sendBit(int type,uint8_t * buffer,int * buffer_len,uint8_t * output_bitstream,int * output_len)
{
	(*buffer) <<= 1;
	(*buffer_len)++;
	if(type==1)		
		(*buffer) |= 0x01;
	else if(type==0)
		(*buffer) |= 0x00;

	if((*buffer_len) == 8)
	{
		output_bitstream[*output_len] = (*buffer);
		(*output_len)++;
		(*buffer) = 0x00;
		(*buffer_len) = 0;
	}

}

void Comp::arithmeticEncoder(int * ori,int ori_len,int range,int * cum_cnt,
							 uint8_t * output_bitstream,int * output_len)
{
	// ori: input sequence
	// range: input sequence range is 0:range
	// ori_len: input sequence length
	// sym_cnt: symbol count
	// buffer: temp buffer
	// buffer_len : temp buffer_len
	// output_bitstream: output bitstream
	// output_len: output length(byte length)

	// initialization
	int * sym_cnt = new int[range+1];
	memset(sym_cnt,0,sizeof(int)*(range+1));

	memset(cum_cnt,0,sizeof(int)*(range+2));
	int total_cnt = ori_len;
	int scale3 = 0;
	// m = 16/24
	uint64_t lower = 0;			  // lower = 0x0000
	uint64_t upper = pow(2,M)-1;   // upper = 0xffff
	uint64_t quarter1 = upper/4 + 1;
	uint64_t middle = quarter1 * 2;
	uint64_t quarter3 = quarter1 * 3;
	uint64_t lower_temp,upper_temp = 0;
	uint8_t buffer = 0x00;
	int buffer_len = 0;

	// calculate cumulate count and total count
	for(int i=0;i<ori_len;i++)
		sym_cnt[ori[i]]++;
	// note:cmu_cnt[0]=0, cum_cnt[i] = num[0]+...+num[i-1]
	for(int i=0;i<range+1;i++)
		for(int j=0;j<=i;j++)
			cum_cnt[i+1] +=sym_cnt[j];

	// start encoding
	for(int i=0;i<ori_len;i++)
	{
		// cout<<i<<" ";
		// update lower and upper
		lower_temp = lower + floor(double((upper-lower+1)*cum_cnt[ori[i]]/total_cnt));
		upper_temp = lower + floor(double((upper-lower+1)*cum_cnt[ori[i]+1]/total_cnt))-1;
		lower = lower_temp;
		upper = upper_temp;

		while((upper<=middle) || (lower>middle) || ((lower>quarter1)&&(upper<=quarter3)))
			//while((upper<=middle) || (lower>middle))
		{
			if(upper<=middle) //E1
			{
				// send MSB and rescale
				Comp::sendBit(0,&buffer,&buffer_len,output_bitstream,output_len);
				lower = lower * 2;
				upper = upper * 2 + 1;
				// consider past E3 condition
				while(scale3>0)
				{
					// send complement of MSB
					Comp::sendBit(1,&buffer,&buffer_len,output_bitstream,output_len);
					scale3--;
				}
			}
			else if(lower>middle) // E2
			{
				// send MSB and rescale
				Comp::sendBit(1,&buffer,&buffer_len,output_bitstream,output_len);
				lower = (lower - middle) * 2;
				upper = (upper - middle) * 2 + 1;
				// consider past E3 condition
				while(scale3>0)
				{
					// send complement of MSB
					Comp::sendBit(0,&buffer,&buffer_len,output_bitstream,output_len);
					scale3--;
				}
			}

			else if(((lower>quarter1)&&(upper<=quarter3)))// E3 condition holds
			{
				// rescale lower and upper
				lower = (lower - quarter1) * 2;
				upper = (upper - quarter1) * 2 + 1;
				scale3++;
			}
		}
	}

	int lower_bin[M];
	// end encoding: send lower and scale3
	lower_temp = lower;
	for(int i=0;i<M;i++)
	{
		lower_bin[M-1-i] = lower_temp%2;
		lower_temp = (lower_temp - lower_temp%2)/2;
	}



	// in order to make sure decoding, mannually add something here:
	for(int i=0;i<20;i++)
	{
		// update lower and upper
		lower_temp = lower + floor(double((upper-lower+1)*cum_cnt[ori[i]]/total_cnt));
		upper_temp = lower + floor(double((upper-lower+1)*cum_cnt[ori[i]+1]/total_cnt))-1;
		lower = lower_temp;
		upper = upper_temp;

		while((upper<=middle) || (lower>middle) || ((lower>quarter1)&&(upper<=quarter3)))
			//while((upper<=middle) || (lower>middle))
		{
			if(upper<=middle) //E1
			{
				// send MSB and rescale
				Comp::sendBit(0,&buffer,&buffer_len,output_bitstream,output_len);
				lower = lower * 2;
				upper = upper * 2 + 1;
				// consider past E3 condition
				while(scale3>0)
				{
					// send complement of MSB
					Comp::sendBit(1,&buffer,&buffer_len,output_bitstream,output_len);
					scale3--;
				}
			}
			else if(lower>middle) // E2
			{
				// send MSB and rescale
				Comp::sendBit(1,&buffer,&buffer_len,output_bitstream,output_len);
				lower = (lower - middle) * 2;
				upper = (upper - middle) * 2 + 1;
				// consider past E3 condition
				while(scale3>0)
				{
					// send complement of MSB
					Comp::sendBit(0,&buffer,&buffer_len,output_bitstream,output_len);
					scale3--;
				}
			}

			else if(((lower>quarter1)&&(upper<=quarter3)))// E3 condition holds
			{
				// rescale lower and upper
				lower = (lower - quarter1) * 2;
				upper = (upper - quarter1) * 2 + 1;
				scale3++;
			}
		}
	}



	sendBit(lower_bin[0],&buffer,&buffer_len,output_bitstream,output_len);
	for(int i=0;i<scale3;i++)
		sendBit(scale3,&buffer,&buffer_len,output_bitstream,output_len);// how large is scale3??
	for(int i=1;i<M;i++)
		sendBit(lower_bin[i],&buffer,&buffer_len,output_bitstream,output_len);



	// end encoding with zeros
	while(buffer_len!=0)
		Comp::sendBit(0,&buffer,&buffer_len,output_bitstream,output_len);


}

void DeComp::getTag(uint8_t * input_bitstream,int * byte_pos,int * bit_pos,uint64_t * tag)
{
	uint8_t this_byte = input_bitstream[(*byte_pos)];
	uint8_t next_byte = input_bitstream[(*byte_pos)+1];
	uint8_t third_byte = input_bitstream[(*byte_pos)+2];
	uint8_t forth_byte = input_bitstream[(*byte_pos)+3];//
	uint8_t fifth_byte = input_bitstream[(*byte_pos)+4];//

	int * tag_bin = new int[M];
	(* tag) = 0;
	int * long_tag_bin = new int[M+8];
	for(int i=0;i<8;i++)
	{
		long_tag_bin[7-i] = this_byte%2;
		long_tag_bin[7-i+8] = next_byte%2;
		long_tag_bin[7-i+16] = third_byte%2;
		long_tag_bin[7-i+24] = forth_byte%2;
		long_tag_bin[7-i+32] = fifth_byte%2;
		this_byte /= 2;
		next_byte /= 2;
		third_byte /= 2;
		forth_byte /= 2;
		fifth_byte /= 2;
		/*
		this_byte = (this_byte - this_byte%2)/2;
		next_byte = (next_byte - next_byte%2)/2;
		third_byte = (third_byte - third_byte%2)/2;
		forth_byte = (forth_byte - forth_byte%2)/2;//
		fifth_byte = (fifth_byte - fifth_byte%2)/2;//
		*/
	}
	for(int i=(*bit_pos);i<(*bit_pos)+M;i++)
	{
		tag_bin[i-(*bit_pos)] = long_tag_bin[i];
		(*tag) = (*tag) + pow(2,M-1-(i-(*bit_pos))) * tag_bin[i-(*bit_pos)];
	}

	delete [] long_tag_bin;
	delete [] tag_bin;
}

void DeComp::arithmeticDecoder(uint8_t * input_bitstream,int input_len,int range,
							   int * cum_cnt,int * ori,int ori_len)
{
	// input_bitstream: input
	// input_len: input length
	// range: original sequence's range: 0 to range
	// cum_cnt: cumulate count
	// ori: output original sequence
	// ori_len: original sequence length(also=total count)
	// m = 16 bits;

	// initialization
	int total_cnt = ori_len;
	uint64_t lower = 0;					// lower	= 0000...0
	uint64_t upper = pow(2,M)-1;		// upper	= 1111...1
	uint64_t quarter1 = upper / 4 + 1;	// quarter1 = 0100...0
	uint64_t middle = quarter1 * 2;		// middle	= 1000...0 
	uint64_t quarter3 = quarter1 * 3;	//quarter3	= 1100...0
	uint64_t lower_temp,upper_temp = 0;
	int * byte_pos = new int;
	int * bit_pos = new int;
	uint64_t * tag = new uint64_t;
	(*byte_pos) = 0;
	(*bit_pos) = 0;
	(*tag) = 0;
	int decode_len = 0;
	int k = 0;
	getTag(input_bitstream,byte_pos,bit_pos,tag);

	// start decoding
	while(decode_len<ori_len)
	{
		// decode a symbol

		k = 0;
		while( floor(double((((*tag)-lower+1)*total_cnt-1)/(upper-lower+1))) >= cum_cnt[k])
			k++;
		ori[decode_len] = k-1;
		decode_len++;
		if(decode_len == ori_len)
			break;

		lower_temp = lower + floor((upper-lower+1)*cum_cnt[ori[decode_len-1]]/total_cnt);//
		upper_temp = lower + floor((upper-lower+1)*cum_cnt[ori[decode_len-1]+1]/total_cnt) - 1;//
		lower = lower_temp;
		upper = upper_temp;

		// rescale lower and upper, update tag
		while((upper<=middle) || (lower>middle) || ((lower>quarter1)&&(upper<=quarter3)))
			//while((upper<=middle) || (lower>middle))
		{
			if(upper<=middle)
			{
				lower = lower * 2;
				upper = upper * 2 + 1;
				(*bit_pos)++;
				if((*bit_pos)==8)
				{
					(*bit_pos) = 0;
					(*byte_pos)++;
				}
				getTag(input_bitstream,byte_pos,bit_pos,tag);
			}
			else if(lower>middle)
			{
				lower = (lower - middle) * 2;
				upper = (upper - middle) * 2 + 1;
				(*bit_pos)++;
				if((*bit_pos)==8)
				{
					(*bit_pos) = 0;
					(*byte_pos)++;
				}
				getTag(input_bitstream,byte_pos,bit_pos,tag);
			}
			else if((lower>quarter1)&&(upper<=quarter3))
			{
				lower = (lower - quarter1) * 2;
				upper = (upper - quarter1) * 2 + 1;
				(*bit_pos)++;
				if((*bit_pos)==8)
				{
					(*bit_pos) = 0;
					(*byte_pos)++;
				}
				getTag(input_bitstream,byte_pos,bit_pos,tag);
				// complement new MSB of tag
				if((*tag)<middle)
					(*tag) = (*tag) + middle;
				else if((*tag)>=middle)
					(*tag) = (*tag) - middle;//???



			}
		}

	}

}




void binChar2hex(unsigned char * binChar, int * hex)
{
	// convert one float
	for(int i=0;i<4;i++)
	{
		hex[i*2+1] = binChar[i]%16;
		hex[i*2] = (binChar[i]-binChar[i]%16)/16;
	}
}

// bin to float: float convertion from hexidecimal to decimal
// binChar: input hexidecimal array
void btof(unsigned char * binChar,float * f)
{
	int bin[32];
	int temp[4];
	float fraction;
	int sign;
	int mantissa;

	temp[0] = (int)binChar[0];
	temp[1] = (int)binChar[1];
	temp[2] = (int)binChar[2];
	temp[3] = (int)binChar[3];
	for(int j=7;j>=0;j--)
	{
		bin[j] = temp[0]%2;
		bin[j+8] = temp[1]%2;
		bin[j+16] = temp[2]%2;
		bin[j+24] = temp[3]%2;
		/*temp[0] = (temp[0] - temp[0]%2)/2;
		temp[1] = (temp[1] - temp[1]%2)/2;
		temp[2] = (temp[2] - temp[2]%2)/2;
		temp[3] = (temp[3] - temp[3]%2)/2;*/

		temp[0] /= 2;
		temp[1] /= 2;
		temp[2] /= 2;
		temp[3] /= 2;
	}
	sign = bin[0];
	mantissa = 0;
	for(int j=1;j<=8;j++)
	{
		mantissa <<= 1;
		mantissa += bin[j];
	}
	fraction = 1;
	for(int j=9;j<=31;j++)
		fraction += bin[j] * pow(0.5,j-8);
	*f = fraction * pow(-1,sign) * pow(2,mantissa-127);
}

// pairsComp
void Comp::pairsComp(FILE * fpW,XMLElement * peaks)
{

	const char * base64char = peaks -> GetText();
	int len = strlen(base64char);
	int flen = 0;
	unsigned char * bin = new unsigned char[len*3/4];
	bin = unbase64(base64char,len,&flen);// credit to ...
	int pairs_len = flen / 4;


	//initialization
	Comp::zero_len = new int [pairs_len/2];
	Comp::diff_reshape = new int [pairs_len * 8];
	Comp::this_mz_hex = new int[8];
	Comp::last_mz_hex = new int[8];
	Comp::diff_reshape_len = 0;

	binChar2hex(bin,Comp::last_mz_hex);

	Comp::pointer = new int [pairs_len];
	Comp::residual = new uint8_t [ 4 * pairs_len];
	Comp::unmatch_data = new uint8_t [4 * (pairs_len)];
	Comp::this_intensity_hex = new int [8];
	Comp::find_intensity_hex = new int [8];
	int match_len = 0;
	int unmatch_len = 0;

	for(int i = 0;i<8;i++)
		Comp::diff_reshape[i] = Comp::last_mz_hex[i];
	Comp::diff_reshape_len = 8;
	Comp::zero_len[0] = 0;

	for(int i=0;i<4;i++)
		Comp::unmatch_data[i] = bin[4+i];
	unmatch_len ++;

	Comp::pointer[0] = 0;



	// decide match type
	// read first 50 intensities, decide match type
	int match_type0_cnt = 0;
	for(int i=2;i<30;i+=2)
	{
		Comp::flag_find_match = false;
		binChar2hex(bin+(i+1)*4,Comp::this_intensity_hex);
		for (int j = i/2-1;j>= (((i/2 - INT_SEARCH_RANGE)<0)?0:(i/2 - INT_SEARCH_RANGE));j--)
		{
			binChar2hex(bin+(j*2+1)*4,Comp::find_intensity_hex);
			bool stoop = true;
			if (Comp::this_intensity_hex[0]==Comp::find_intensity_hex[0] &&
				Comp::this_intensity_hex[1]==Comp::find_intensity_hex[1] &&
				Comp::this_intensity_hex[4]==Comp::find_intensity_hex[4] &&
				Comp::this_intensity_hex[5]==Comp::find_intensity_hex[5] &&
				Comp::this_intensity_hex[6]==Comp::find_intensity_hex[6] &&
				Comp::this_intensity_hex[7]==Comp::find_intensity_hex[7])
			{
				match_type0_cnt++;
				break;
			}
		}

		if(match_type0_cnt>5)
			Comp::match_type = 0;
		else
			Comp::match_type = 1;
	}




	// start encoding

	if(Comp::match_type == 0)
	{
		// type 0: match 0:1,4:7
		for(int i = 2;i < pairs_len;i += 2)
		{
			// mz compression
			binChar2hex(bin+i*4,Comp::this_mz_hex);
			Comp::cnt = 0;
			Comp::flag_cnt_zero = true;

			for(int j = 0;j < 8;j++)
			{
				Comp::diff[j] = (Comp::this_mz_hex[j] - Comp::last_mz_hex[j] + 16)%16;
				if (Comp::diff[j]!=0)
					Comp::flag_cnt_zero = false;
				if (Comp::flag_cnt_zero)
					Comp::cnt++;
				else
				{				
					Comp::diff_reshape[Comp::diff_reshape_len] = Comp::diff[j];
					Comp::diff_reshape_len++;
				}
			}
			Comp::zero_len[i/2] = Comp::cnt;
			binChar2hex(bin+i*4,Comp::last_mz_hex);


			// intensity compression
			Comp::flag_find_match = false;
			binChar2hex(bin+(i+1)*4,Comp::this_intensity_hex);

			for (int j = i/2-1;j>= (((i/2 - INT_SEARCH_RANGE)<0)?0:(i/2 - INT_SEARCH_RANGE));j--)
			{
				binChar2hex(bin+(j*2+1)*4,Comp::find_intensity_hex);
				bool stoop = true;

				if (Comp::this_intensity_hex[0]==Comp::find_intensity_hex[0] &&
					Comp::this_intensity_hex[1]==Comp::find_intensity_hex[1] &&
					Comp::this_intensity_hex[4]==Comp::find_intensity_hex[4] &&
					Comp::this_intensity_hex[5]==Comp::find_intensity_hex[5] &&
					Comp::this_intensity_hex[6]==Comp::find_intensity_hex[6] &&
					Comp::this_intensity_hex[7]==Comp::find_intensity_hex[7])
				{
					Comp::flag_find_match = true;
					Comp::pointer[i/2] = i/2 - j;
					Comp::residual[match_len] = (Comp::this_intensity_hex[2]<<4) + Comp::this_intensity_hex[3];
					match_len++;
					break;
				}
			}
			if (Comp::flag_find_match == false)
			{
				Comp::pointer[i/2] = 0;
				for(int i_unmatch=0;i_unmatch<4;i_unmatch++)
					Comp::unmatch_data[unmatch_len*4+i_unmatch] = bin[4*(i+1)+i_unmatch];
				unmatch_len++;
			}
		}
	}

	else if(Comp::match_type == 1)
	{
		// type 1: match 0:1
		for(int i = 2;i < pairs_len;i += 2)
		{
			// mz compression
			binChar2hex(bin+i*4,Comp::this_mz_hex);
			Comp::cnt = 0;
			Comp::flag_cnt_zero = true;

			for(int j = 0;j < 8;j++)
			{
				Comp::diff[j] = (Comp::this_mz_hex[j] - Comp::last_mz_hex[j] + 16)%16;
				if (Comp::diff[j]!=0)
					Comp::flag_cnt_zero = false;
				if (Comp::flag_cnt_zero)
					Comp::cnt++;
				else
				{				
					Comp::diff_reshape[Comp::diff_reshape_len] = Comp::diff[j];
					Comp::diff_reshape_len++;
				}
			}
			Comp::zero_len[i/2] = Comp::cnt;
			binChar2hex(bin+i*4,Comp::last_mz_hex);


			// intensity compression
			Comp::flag_find_match = false;
			binChar2hex(bin+(i+1)*4,Comp::this_intensity_hex);

			for (int j = i/2-1;j>= (((i/2 - INT_SEARCH_RANGE)<0)?0:(i/2 - INT_SEARCH_RANGE));j--)
			{
				binChar2hex(bin+(j*2+1)*4,Comp::find_intensity_hex);
				bool stoop = true;

				if (Comp::this_intensity_hex[0]==Comp::find_intensity_hex[0] &&
					Comp::this_intensity_hex[1]==Comp::find_intensity_hex[1])
				{
					Comp::flag_find_match = true;
					Comp::pointer[i/2] = i/2 - j;
					Comp::residual[match_len] = (Comp::this_intensity_hex[2]<<4) + Comp::this_intensity_hex[3];
					Comp::residual[match_len+1] = (Comp::this_intensity_hex[4]<<4) + Comp::this_intensity_hex[5];
					Comp::residual[match_len+2] = (Comp::this_intensity_hex[6]<<4) + Comp::this_intensity_hex[7];
					match_len+=3;
					break;
				}
			}
			if (Comp::flag_find_match == false)
			{
				Comp::pointer[i/2] = 0;
				for(int i_unmatch=0;i_unmatch<4;i_unmatch++)
					Comp::unmatch_data[unmatch_len*4+i_unmatch] = bin[4*(i+1)+i_unmatch];
				unmatch_len++;
			}
		}
	}








	// in case of odd length of diff_reshape
	Comp::diff_reshape[Comp::diff_reshape_len] = 0;
	Comp::diff_reshape[Comp::diff_reshape_len+1] = 0;


	// save 2 diff_reshape[i] in 1 byte
	Comp::diff_reshape_byte = new uint8_t[(int)ceil((double)Comp::diff_reshape_len/2)];
	for(int i=0;i<(int)ceil((double)Comp::diff_reshape_len/2);i++)
		Comp::diff_reshape_byte[i] = (Comp::diff_reshape[i*2]<<4) + Comp::diff_reshape[i*2+1];


	// block 1: pairs_len
	fwrite(&pairs_len,sizeof(int),1,fpW);

	// block 2: m/z 
	// arithmetic coding for zero_len
	Comp::output_bitstream = new uint8_t[pairs_len];
	Comp::output_len = new int;
	*(Comp::output_len) = 0;
	int arith_range = 8; //0:8
	Comp::cum_cnt = new int[arith_range+2];

	fwrite(&arith_range,sizeof(int),1,fpW);
	Comp::arithmeticEncoder(Comp::zero_len,pairs_len/2,arith_range,
		Comp::cum_cnt,Comp::output_bitstream,Comp::output_len);
	fwrite(Comp::cum_cnt,sizeof(int),arith_range+2,fpW);
	fwrite(Comp::output_len,sizeof(int),1,fpW);
	fwrite(Comp::output_bitstream,sizeof(uint8_t),(*Comp::output_len),fpW);

	// diff_reshape
	fwrite(&Comp::diff_reshape_len,sizeof(int),1,fpW);
	fwrite(Comp::diff_reshape_byte,sizeof(uint8_t),(int)ceil((double)Comp::diff_reshape_len/2),fpW);


	// check error
	/*std::cout<<"check error"<<std::endl;

	std::cout<<std::endl<<std::endl<<"zero_len"<<std::endl;
	for(int i=0;i<pairs_len/2;i++)
	std::cout<<Comp::zero_len[i]<<" ";

	std::cout<<std::endl<<std::endl<<"diff_reshape"<<std::endl;
	for(int i=0;i<Comp::diff_reshape_len;i++)
	std::cout<<Comp::diff_reshape[i]<<" ";*/



	// block 3: intensity
	// match_type
	fwrite(&Comp::match_type,sizeof(int),1,fpW);
	// arithmetic coding for pointer
	Comp::output_bitstream = new uint8_t[pairs_len];
	*(Comp::output_len) = 0;
	arith_range = INT_SEARCH_RANGE; //0:INT_SEARCH_RANGE

	fwrite(&arith_range,sizeof(int),1,fpW);
	Comp::cum_cnt = new int[arith_range+2];
	Comp::arithmeticEncoder(Comp::pointer,pairs_len/2,arith_range,
		Comp::cum_cnt,Comp::output_bitstream,Comp::output_len);
	fwrite(Comp::cum_cnt,sizeof(int),arith_range+2,fpW);
	fwrite(Comp::output_len,sizeof(int),1,fpW);
	fwrite(Comp::output_bitstream,sizeof(uint8_t),(*Comp::output_len),fpW);

	// residual
	fwrite(&match_len,sizeof(int),1,fpW);
	fwrite(Comp::residual,sizeof(uint8_t),match_len,fpW);

	// unmatch_data
	fwrite(&unmatch_len,sizeof(int),1,fpW);
	fwrite(Comp::unmatch_data,sizeof(uint8_t),unmatch_len*4,fpW);



	// check error
	/*std::cout<<std::endl<<std::endl<<"pointer"<<std::endl;
	for(int i=0;i<pairs_len/2;i++)
	std::cout<<Comp::pointer[i]<<" ";
	std::cout<<std::endl<<std::endl<<"residual"<<std::endl;
	for(int i=0;i<match_len;i++)
	std::cout<<(int)Comp::residual[i]<<" ";
	std::cout<<std::endl<<std::endl<<"unmatch_data"<<std::endl;
	for(int i=0;i<unmatch_len*4;i++)
	std::cout<<(int)Comp::unmatch_data[i]<<" ";*/


	// free the space
	delete [] Comp::this_mz_hex;
	delete [] Comp::last_mz_hex;
	delete [] Comp::zero_len;
	delete [] Comp::diff_reshape;
	delete [] Comp::diff_reshape_byte;

	delete [] Comp::this_intensity_hex;
	delete [] Comp::find_intensity_hex;
	delete [] Comp::pointer;
	delete [] Comp::residual;
	delete [] Comp:: unmatch_data;

	delete [] Comp:: output_bitstream;
	delete [] Comp:: output_len;

	delete [] bin;
	delete [] Comp::cum_cnt;


}

void DeComp::pairsDecomp(FILE ** fp, XMLElement * scan, XMLDocument * doc)
{

	fread(&pairs_len,sizeof(int),1,*fp);

	if(pairs_len<50)
		return;

	// decompress mz block
	// read data
	fread(&range,sizeof(int),1,*fp);
	cum_cnt = new int[range + 2];
	fread(cum_cnt,sizeof(int),range+2,*fp);
	fread(&input_len,sizeof(int),1,*fp);
	input_bitstream = new uint8_t[input_len];
	fread(input_bitstream,sizeof(uint8_t),input_len,*fp);

	zero_len = new int[pairs_len/2];
	DeComp::arithmeticDecoder(input_bitstream,input_len,
		range,cum_cnt,zero_len,pairs_len/2);
	fread(&diff_reshape_len,sizeof(int),1,*fp);
	uint8_t * diff_reshape = new uint8_t[(int)ceil((double)diff_reshape_len/2)];
	fread(diff_reshape,sizeof(uint8_t),(int)ceil((double)diff_reshape_len/2),*fp);

	diff_reshape_hex = new int[diff_reshape_len+10];
	for(int i=0;i<(int)ceil((double)diff_reshape_len/2);i++)
	{
		diff_reshape_hex[i*2] = (diff_reshape[i] - diff_reshape[i]%16)/16;
		diff_reshape_hex[i*2+1] = diff_reshape[i]%16;
	}



	// check error
	/*std::cout<<"check error"<<std::endl;

	std::cout<<std::endl<<std::endl<<"zero_len"<<std::endl;
	for(int i=0;i<pairs_len/2;i++)
	std::cout<<zero_len[i]<<" ";

	std::cout<<std::endl<<std::endl<<"diff_reshape"<<std::endl;
	for(int i=0;i<diff_reshape_len;i++)
	std::cout<<diff_reshape_hex[i]<<" ";*/



	last_mz_hex = new int[8];
	this_mz_hex = new int[8];
	diff_hex = new int[8];
	bin = new unsigned char[pairs_len*4];
	pos = 0;



	for(int i=0;i<8;i++)
		last_mz_hex[i] = diff_reshape_hex[i];
	bin[0] = last_mz_hex[0]*16 + last_mz_hex[1];
	bin[1] = last_mz_hex[2]*16 + last_mz_hex[3];
	bin[2] = last_mz_hex[4]*16 + last_mz_hex[5];
	bin[3] = last_mz_hex[6]*16 + last_mz_hex[7];
	pos += 8;

	/*
	// check error
	std::cout<<std::endl<<"mz in hex"<<std::endl;
	*/
	// conjacent data
	for(int i=1;i<pairs_len/2;i++)
	{
		// decompress this_mz_hex
		for(int j=0;j<zero_len[i];j++)
		{
			diff_hex[j] = 0;
			this_mz_hex[j] = (diff_hex[j] + last_mz_hex[j])%16;
		}
		for(int j=zero_len[i];j<8;j++)
		{
			diff_hex[j] = diff_reshape_hex[pos+j-zero_len[i]];
			this_mz_hex[j] = (diff_hex[j] + last_mz_hex[j])%16;
		}
		pos += (8 - zero_len[i]);

		// put this_mz_hex in bin
		for(int j=0;j<4;j++)
			bin[8*i+j] = this_mz_hex[j*2]*16 + this_mz_hex[j*2+1];

		for(int j=0;j<8;j++)
			last_mz_hex[j] = this_mz_hex[j];


		// check error
		//for(int cout=0;cout<8;cout++)
		//std::cout<<this_mz_hex[cout]<<" ";
	}



	// decompress intensity block
	// read data
	fread(&match_type,sizeof(int),1,*fp);
	// pointer
	fread(&range,sizeof(int),1,*fp);
	cum_cnt = new int[range + 2];
	fread(cum_cnt,sizeof(int),range+2,*fp);
	fread(&input_len2,sizeof(int),1,*fp);
	input_bitstream2 = new uint8_t[input_len2];
	fread(input_bitstream2,sizeof(uint8_t),input_len2,*fp);

	pointer = new int[pairs_len/2];
	//int * zero_len = new int[pairs_len/2];
	DeComp::arithmeticDecoder(input_bitstream2,input_len,
		range,cum_cnt,pointer,pairs_len/2);

	// residual
	fread(&match_len,sizeof(int),1,*fp);
	residual = new uint8_t[match_len];
	fread(residual,sizeof(uint8_t),match_len,*fp);

	// unmatch_data
	fread(&unmatch_len,sizeof(int),1,*fp);
	unmatch_data = new uint8_t[unmatch_len*4];
	fread(unmatch_data,sizeof(uint8_t),unmatch_len*4,*fp);


	// check error
	/*std::cout<<std::endl<<std::endl<<"pointer"<<std::endl;
	for(int i=0;i<pairs_len/2;i++)
	std::cout<<pointer[i]<<" ";
	std::cout<<std::endl<<std::endl<<"residual"<<std::endl;
	for(int i=0;i<match_len;i++)
	std::cout<<(int)residual[i]<<" ";
	std::cout<<std::endl<<std::endl<<"unmatch_data"<<std::endl;
	for(int i=0;i<unmatch_len*4;i++)
	std::cout<<(int)unmatch_data[i]<<" ";*/



	unmatch_cnt = 0;
	match_cnt = 0;
	find_int_hex = new int[8];
	this_int_hex = new int[8];
	// conjacent data

	if(match_type == 0)
	{
		// type 0: match 0:1,4:7

		for(int i=0;i<pairs_len/2;i++)
		{
			if(pointer[i]==0)
			{
				for(int j=0;j<4;j++)
					bin[i*8+j+4] = unmatch_data[unmatch_cnt*4+j];
				unmatch_cnt++;
			}
			else
			{
				this_int_hex[3] = residual[match_cnt]%16;
				this_int_hex[2] = (residual[match_cnt] - residual[match_cnt]%16)/16;
				for(int j=0;j<4;j++)
				{
					find_int_hex[j*2+1] = bin[8*(i-pointer[i])+4+j]%16;
					find_int_hex[j*2] = (bin[8*(i-pointer[i])+4+j] - bin[8*(i-pointer[i])+4+j]%16)/16;
				}
				this_int_hex[0] = find_int_hex[0];
				this_int_hex[1] = find_int_hex[1];
				this_int_hex[4] = find_int_hex[4];
				this_int_hex[5] = find_int_hex[5];
				this_int_hex[6] = find_int_hex[6];
				this_int_hex[7] = find_int_hex[7];
				match_cnt++;

				for(int j=0;j<4;j++)
					bin[i*8+j+4] = this_int_hex[2*j]*16 + this_int_hex[2*j+1];

			}
		}

	}
	else if(match_type == 1)
	{
		// type 1: match 0:1
		for(int i=0;i<pairs_len/2;i++)
		{
			if(pointer[i]==0)
			{
				for(int j=0;j<4;j++)
					bin[i*8+j+4] = unmatch_data[unmatch_cnt*4+j];
				unmatch_cnt++;
			}
			else
			{
				this_int_hex[3] = residual[match_cnt]%16;
				this_int_hex[2] = (residual[match_cnt] - residual[match_cnt]%16)/16;
				this_int_hex[5] = residual[match_cnt+1]%16;
				this_int_hex[4] = (residual[match_cnt+1] - residual[match_cnt+1]%16)/16;
				this_int_hex[7] = residual[match_cnt+2]%16;
				this_int_hex[6] = (residual[match_cnt+2] - residual[match_cnt+2]%16)/16;

				for(int j=0;j<4;j++)
				{
					find_int_hex[j*2+1] = bin[8*(i-pointer[i])+4+j]%16;
					find_int_hex[j*2] = (bin[8*(i-pointer[i])+4+j] - bin[8*(i-pointer[i])+4+j]%16)/16;
				}
				this_int_hex[0] = find_int_hex[0];
				this_int_hex[1] = find_int_hex[1];
				match_cnt+=3;

				for(int j=0;j<4;j++)
					bin[i*8+j+4] = this_int_hex[2*j]*16 + this_int_hex[2*j+1];

			}
		}

	}








	// check error
	/*std::cout<<std::endl<<"bin for decompressed:"<<std::endl;
	for(int i=0;i<pairs_len*4;i++)
	std::cout<<(int)bin[i]<<" ";
	std::cout<<std::endl;*/






	// base64 code for pairs
	// bin to base64char
	flen = 0;
	base64char = base64( bin,4*pairs_len,&flen);


	// add child "peaks" to element "scan"
	XMLElement * peaks = scan -> FirstChildElement("peaks");
	peaks -> SetText(base64char);
	//peaks -> SetText(base64char);
	//peaks -> SetAttribute("precision",32);
	//peaks -> SetAttribute("byteOrder","network");
	//peaks -> SetAttribute("pairOrder","m/z-int");
	//scan -> LinkEndChild(peaks);
	// attributes
	// text


	delete [] bin;
	delete [] base64char;

	delete [] cum_cnt;
	delete [] input_bitstream;
	delete [] input_bitstream2;
	delete [] diff_reshape;
	delete [] residual;
	delete [] unmatch_data;
	delete [] zero_len;
	delete [] pointer;

	delete [] this_mz_hex;
	delete [] last_mz_hex;
	delete [] diff_hex;
	delete [] this_int_hex;
	delete [] find_int_hex;
	delete [] diff_reshape_hex;
}

void FPMSComp(std::string input_path, std::string output_folder)
{
	std::size_t pos_last = input_path.find_last_of("\\");
	std::size_t pos_ext = input_path.find_last_of(".mzXML");
	std::string filename = input_path.substr(pos_last+1,pos_ext-6-pos_last);
	std::string output_pairs = "." + output_folder + "\\pairsCompressed_" + filename + ".bin";
	std::string output_xml = "." + output_folder + "\\struct_" + filename + ".xml";



	//std::cout<<"start compressing"<<std::endl<<input_path<<std::endl;

	XMLDocument doc;
	XMLError eResult = doc.LoadFile(input_path.c_str());
	XMLElement * mzXML = doc.RootElement();
	XMLElement * msRun = mzXML -> FirstChildElement("msRun");
	XMLElement * scan = msRun -> FirstChildElement("scan");
	XMLElement * peaks = scan -> FirstChildElement("peaks");
	int scanCount = 0;
	msRun->QueryIntAttribute("scanCount",&scanCount);

	FILE * fpW = fopen(output_pairs.c_str(),"wb");
	fwrite(&scanCount,sizeof(int),1,fpW);
	int peaksCount;

	for(int cnt=0;cnt<scanCount;cnt++)
	{
		peaks = scan -> FirstChildElement("peaks");
		// pairsComp and delete
		scan -> QueryIntAttribute("peaksCount",&peaksCount);
		if(peaksCount>=50)
		{
			Comp::pairsComp(fpW,peaks);
			peaks -> SetText("");
		}
		//std::cout<<cnt+1<<" ";
		if(cnt<scanCount-1)
			scan = scan -> NextSiblingElement();
	}

	eResult = doc.SaveFile(output_xml.c_str());


	fclose((fpW));

	//std::cout<<"end compressing"<<std::endl<<std::endl;

}

void FPMSDecomp(std::string input_path ,std::string output_folder)
{

	/*char output_path[200];
	strcpy(output_path,folder_out);
	strcat(output_path,"PairsCompressed_");
	strcat(output_path,filename);

	char output_path_struct[200];
	strcpy(output_path_struct,folder_out);
	strcat(output_path_struct,"FileStruct_");
	strcat(output_path_struct,filename);
	strcat(output_path_struct,".xml");


	char decompress_path[200];
	strcpy(decompress_path,folder_out);
	strcat(decompress_path,"Decompressed_");
	strcat(decompress_path,filename);
	strcat(decompress_path,".mzxml");
	*/


	std::string path_pairs = input_path;

	std::size_t pos_last = input_path.find_last_of("\\");
	std::size_t pos_ext = input_path.find_last_of(".bin");
	std::string filename = input_path.substr(pos_last+17,pos_ext-4-pos_last-16);
	std::string compressed_folder = input_path.substr(0,pos_last);
	std::string path_struct = compressed_folder + "\\" + "Struct_" + filename + ".xml";

	std::string path_decompress = "." + output_folder + "\\" + "Decompressed_" + filename + ".mzXML";

	XMLDocument doc;
	doc.LoadFile(path_struct.c_str());
	XMLElement * mzXML = doc.RootElement();
	XMLElement * msRun = mzXML -> FirstChildElement("msRun");
	XMLElement * scan = msRun -> FirstChildElement("scan");
	//XMLElement * peaks = scan -> FirstChildElement("peaks");
	int scanCount;
	int peaksCount;

	FILE * fpR;
	fpR = fopen(path_pairs.c_str(),"rb");
	fread(&scanCount,sizeof(int),1,fpR);

	std::cout<<std::endl<<"start decompressing"<<std::endl;
	for(int i=0;i<scanCount;i++)
	{
		scan ->QueryIntAttribute("peaksCount",&peaksCount);
		if(peaksCount<50)
		{}
		else
			DeComp::pairsDecomp(&fpR,scan,&doc);
		scan = scan->NextSiblingElement("scan");
		//std::cout<<i<<" ";
	}

	//std::cout<<std::endl;
	//std::cout<<"end decompressing"<<std::endl;

	doc.SaveFile(path_decompress.c_str());
	fclose(fpR);
}


void CompCmp(char * folder_in,char * folder_out,char * filename)
{
	char ori_filename[200];
	strcpy(ori_filename,folder_in);
	strcat(ori_filename,filename);
	strcat(ori_filename,".mzXML");

	char dc_filename[200];
	strcpy(dc_filename,folder_out);
	strcat(dc_filename,"Decompressed_");
	strcat(dc_filename,filename);
	strcat(dc_filename,".mzXML");


	std::cout<<std::endl<<"start comparing"<<std::endl;

	XMLDocument ori_doc;
	XMLError eResult_ori = ori_doc.LoadFile(ori_filename);
	XMLElement * msRun_ori = ori_doc.RootElement()->FirstChildElement("msRun");
	XMLElement * scan_ori = msRun_ori -> FirstChildElement("scan");
	XMLElement * peaks_ori;
	int scanCount_ori = 0;
	msRun_ori->QueryIntAttribute("scanCount",&scanCount_ori);

	XMLDocument dc_doc;
	XMLError eResult_dc = dc_doc.LoadFile(dc_filename);
	XMLElement * scan_dc = dc_doc.RootElement() -> FirstChildElement("msRun") ->
		FirstChildElement("scan");
	XMLElement * peaks_dc;
	int peaksCount_ori;



	for(int i=0;i<scanCount_ori;i++)
	{
		peaks_ori = scan_ori -> FirstChildElement("peaks");
		peaks_dc = scan_dc -> FirstChildElement("peaks");
		scan_ori -> QueryIntAttribute("peaksCount",&peaksCount_ori);
		if(peaksCount_ori==0)
			continue;

		const char * base64char_ori = peaks_ori -> GetText();
		int len_ori = strlen(base64char_ori);
		int flen_ori = 0;
		unsigned char * bin_ori = new unsigned char[len_ori*3/4];
		bin_ori = unbase64(base64char_ori,len_ori,&flen_ori);

		const char * base64char_dc = peaks_dc -> GetText();
		int len_dc = strlen(base64char_dc);
		int flen_dc = 0;
		unsigned char * bin_dc = new unsigned char[len_dc*3/4];
		bin_dc = unbase64(base64char_dc,len_dc,&flen_dc);

		// cmp
		for(int j=0;j<len_ori*3/4;j++)
		{
			if(bin_ori[j] == bin_dc[j])
				continue;
			else
			{
				std::cout<<"wrong decompression at 'scanNum="<<i+1<<"', 'binNum="<<j+1<<"'"
					<<std::endl<<"original bin is "<<(int)bin_ori[j]<<" and decompressed bin is"<<(int)bin_dc[j]
				<<std::endl<<std::endl;
				break;
			}
		}
		// turn to next scan
		scan_ori = scan_ori -> NextSiblingElement("scan");
		scan_dc = scan_dc -> NextSiblingElement("scan");
	}

	std::cout<<"end comparing"<<std::endl;
}


int main()
{ 	

	/*	comp	*/

	system("dir/s/b *.mzXML > mzXML_dir.txt");
	// input the MassIVE id
	std::string input_folder_path;
	//std::string var[3] = {
	//"\\input\\MSV000080896\\peak\\Data_mzXML",
	//"\\input\\MSV000080905\\mzXML\\Plate1",
	//"\\input\\MSV000080905\\mzXML\\Plate2"};
	//"\\input\\MSV000080905\\mzXML\\Plate3",
	//"\\input\\MSV000080905\\mzXML\\Plate4",
	//"\\input\\MSV000081123\\peaks"};
	// \input\MSV000080905\mzXML\Plate2
	// \input\MSV000080905\mzXML\Plate3
	// \input\MSV000080905\mzXML\Plate4
	// \input\MSV000081123\peaks};
	std::cout<<"please input the file path: "<<std::endl;
	std::cin>>input_folder_path;


	//clock_t start,finish;
	//double totaltime;





	//start=clock();

	std::ifstream fpdir;
	fpdir.open("mzXML_dir.txt");
	std::string file_path;

	// create output folder
	char proj_path[MAX_PATH];   
	getcwd(proj_path,MAX_PATH);
	//std::string output_folder = proj_path;
	std::string output_folder = "\\output" + input_folder_path;
	std::string create_path_cmd = "md ." + output_folder;
	//std::string create_path_cmd = "mkdir -p" + output_folder;
	system(create_path_cmd.c_str());


	// find mzXML files corresponding with MassIVE id
	for(int i=0;!fpdir.eof();i++)
	{
		getline(fpdir,file_path);
		std::size_t pos = file_path.find(input_folder_path);
		//std::size_t pos = file_path.find("Plate1");
		//std::size_t pos = file_path.find("Plate2");
		//std::size_t pos = file_path.find("Plate3");
		//std::size_t pos = file_path.find("Plate4");

		if(pos == std::string::npos)// not this massIVE id
			continue;
		FPMSComp(file_path,output_folder);
	}
	fpdir.close();

	//finish=clock();
	//totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	//std::cout<<"run time"<<totaltime<<"seconds"<<std::endl;





	/*	decomp	*/

	system("dir/s/b *.bin > decomp_dir.txt");
	// input the MassIVE id
	std::string de_input_folder_path;
	//std::string var[6] = {
	//"\\output\\input\\MSV000080896\\peak\\Data_mzXML",
	//"\\output\\input\\MSV000080905\\mzXML\\Plate1",
	//"\\output\\input\\MSV000080905\\mzXML\\Plate2",
	//"\\output\\input\\MSV000080905\\mzXML\\Plate3",
	//"\\output\\input\\MSV000080905\\mzXML\\Plate4",
	//"\\output\\input\\MSV000081123\\peaks"};



	//clock_t start,finish;
	//double totaltime;





	//start=clock();

	std::ifstream de_fpdir;
	de_fpdir.open("decomp_dir.txt");
	//std::string file_path;

	// create output folder
	char de_proj_path[MAX_PATH];   
	getcwd(de_proj_path,MAX_PATH);
	//std::string output_folder = proj_path;
	std::string de_output_folder = "\\decompress" + de_input_folder_path;
	std::string de_create_path_cmd = "md ." + de_output_folder;
	//std::string create_path_cmd = "mkdir -p" + output_folder;
	system(de_create_path_cmd.c_str());


	// find mzXML files corresponding with MassIVE id
	for(int i=0;!fpdir.eof();i++)
	{
		getline(fpdir,file_path);
		std::size_t pos = file_path.find(de_input_folder_path);
		if(pos == std::string::npos)// not this massIVE id
			continue;
		FPMSDecomp(file_path,de_output_folder);
	}
	de_fpdir.close();

	//finish=clock();
	//totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	//std::cout<<"run time"<<totaltime<<"seconds"<<std::endl;




	return 0;
}