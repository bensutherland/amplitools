## Troubleshooting Support ##
Use the following at your own risk, but the solutions should be searchable online and point you in the right direction.     

### Too many open files?   
When running samtools merge with a bamlist, if you have many individual samples being merged together, you may come across the following error:    
```
[E::hts_open_format] Failed to open file "13_mapped_mhap/IonCode_1320_rawlib.sorted.bam" : Too many open files
```

This issue is occurring because there is a soft limit on the number of files that you can have open at one time on your computer. There is also a hard limit. The soft limit can be adjusted, but within reason. Look into these details for your own computer and OS, but Google suggests for a macOS that the following should suffice and be OK for the computer:    

`ulimit -n 10240`       
 
If you have more than 10,000 files (nice!), it will likely be required to look at the specific details of your computer.    

Note: this will reset on a new shell spawn.    

