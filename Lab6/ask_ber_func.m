function [ber,numBits]=ask_ber_func(EbNo, maxNumErrs, maxNumBits)
import com.mathworks.toolbox.comm.BERTool;
totErr=0;
numBits=0;

k=4;
Nsymb=2000;
nsamp=512;

while((totErr<maxNumErrs) && (numBits<maxNumBits))
    if(BERTool.getSimulationStop)
        break;
    end

errors=fsk_errors_noncoh(k,Nsymb,nsamp,EbNo);
totErr=totErr+errors;
numBits=numBits+k*Nsymb;
end
ber=totErr/numBits;