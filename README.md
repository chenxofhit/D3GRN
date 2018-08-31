# D3GRN

D3GRN is a  Data Driven Dynamic network construction method to infer Gene Regulatory Networks.

## System Environment
Ubuntu/Linux Bash, Matlab 2015b.

## Running Step

Step 1,  DREAM4  Network1-5:

```bash
nohup  bash ./main_arni_dream4_multifactorial_network1_5.sh  &
```

the results of bootstraping runs will be generated for  further evaluation.   

Or if you want to run the algorithm for  DREAM5 Network 1 (sudo privilege perhaps needed according to your running environment):

```{bash}
nohup bash ./main_arni_dream5_network1.sh &
```

As the DREAM5 networks are large in scale so it will be time consuming,  so the task is  recommended to run in the background mode. 

Step 2,  Evaluation:

It is recommended to run the following Matlab script in  your desktop IDE, the command will calculate the AUPR as well as plot the curves.

- For  DREAM4  Network1-5:

```bash
infer_grn_arni_dream4_multifactorial_para.m
```

- For DREAM5 Network 1:

```bash
infer_grn_arni_dream5_para_Network1.m
```

Any question, please do not hesitate to  contact me with following address with bash command for decryption:

```{bash} 
echo "Y2hlbnhvZmhpdEBnbWFpbC5jb20K"|base64 -d (use -D instead for Mac users)
```
or [submit issue](https://github.com/chenxofhit/d3grn/issues) in the repository directly.