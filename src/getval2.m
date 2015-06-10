function val=getval2(keys,vals,keyin,itype,nv)
val=[];
for i=1:nv
   k1 = strfind(keys{i},keyin);
   if numel(k1) > 0
      k1 = k1(1);
   else
      k1 = 0;
   end
   if k1 > 0
      str=vals{i};
      k2 = strfind(str,'%');
      if numel(k2) > 0
         k2 = k2(1)-1;
      else
         k2 = numel(str);
      end
      switch itype
         case 's'
            val = strtrim(str(1:k2));
            fprintf(1,'%20s = %s\n',keyin,val);
         case 'f'
            val = str2double(str(1:k2));
            fprintf(1,'%20s = %f\n',keyin,val);
         case 'u'
            val1 = str2double(str(1:k2));
            val = cast(val1,'uint8');
            fprintf(1,'%20s = %x\n',keyin,val);
         case 'b'
            val = cast(bin2dec(str(1:k2)),'uint8');
            fprintf(1,'%20s = %8s\n',keyin,dec2bin(val,8));
         otherwise
            error('unknown itype');
      end
   end
end
return
