glabel func_8099E784
/* 00EB4 8099E784 27BDFF30 */  addiu   $sp, $sp, 0xFF30           ## $sp = FFFFFF30
/* 00EB8 8099E788 AFBF0064 */  sw      $ra, 0x0064($sp)           
/* 00EBC 8099E78C AFBE0060 */  sw      $s8, 0x0060($sp)           
/* 00EC0 8099E790 AFB7005C */  sw      $s7, 0x005C($sp)           
/* 00EC4 8099E794 AFB60058 */  sw      $s6, 0x0058($sp)           
/* 00EC8 8099E798 AFB50054 */  sw      $s5, 0x0054($sp)           
/* 00ECC 8099E79C AFB40050 */  sw      $s4, 0x0050($sp)           
/* 00ED0 8099E7A0 AFB3004C */  sw      $s3, 0x004C($sp)           
/* 00ED4 8099E7A4 AFB20048 */  sw      $s2, 0x0048($sp)           
/* 00ED8 8099E7A8 AFB10044 */  sw      $s1, 0x0044($sp)           
/* 00EDC 8099E7AC AFB00040 */  sw      $s0, 0x0040($sp)           
/* 00EE0 8099E7B0 F7BC0038 */  sdc1    $f28, 0x0038($sp)          
/* 00EE4 8099E7B4 F7BA0030 */  sdc1    $f26, 0x0030($sp)          
/* 00EE8 8099E7B8 F7B80028 */  sdc1    $f24, 0x0028($sp)          
/* 00EEC 8099E7BC F7B60020 */  sdc1    $f22, 0x0020($sp)          
/* 00EF0 8099E7C0 F7B40018 */  sdc1    $f20, 0x0018($sp)          
/* 00EF4 8099E7C4 8CB20000 */  lw      $s2, 0x0000($a1)           ## 00000000
/* 00EF8 8099E7C8 8CAE1C44 */  lw      $t6, 0x1C44($a1)           ## 00001C44
/* 00EFC 8099E7CC 0080B825 */  or      $s7, $a0, $zero            ## $s7 = 00000000
/* 00F00 8099E7D0 00A0F025 */  or      $s8, $a1, $zero            ## $s8 = 00000000
/* 00F04 8099E7D4 3C06809A */  lui     $a2, %hi(D_8099EBB0)       ## $a2 = 809A0000
/* 00F08 8099E7D8 24C6EBB0 */  addiu   $a2, $a2, %lo(D_8099EBB0)  ## $a2 = 8099EBB0
/* 00F0C 8099E7DC 27A4009C */  addiu   $a0, $sp, 0x009C           ## $a0 = FFFFFFCC
/* 00F10 8099E7E0 240701D8 */  addiu   $a3, $zero, 0x01D8         ## $a3 = 000001D8
/* 00F14 8099E7E4 02402825 */  or      $a1, $s2, $zero            ## $a1 = 00000000
/* 00F18 8099E7E8 0C031AB1 */  jal     Graph_OpenDisp              
/* 00F1C 8099E7EC AFAE00B0 */  sw      $t6, 0x00B0($sp)           
/* 00F20 8099E7F0 0C024F46 */  jal     func_80093D18              
/* 00F24 8099E7F4 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 00F28 8099E7F8 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 00F2C 8099E7FC 3C18E700 */  lui     $t8, 0xE700                ## $t8 = E7000000
/* 00F30 8099E800 3C08FA00 */  lui     $t0, 0xFA00                ## $t0 = FA000000
/* 00F34 8099E804 244F0008 */  addiu   $t7, $v0, 0x0008           ## $t7 = 00000008
/* 00F38 8099E808 AE4F02D0 */  sw      $t7, 0x02D0($s2)           ## 000002D0
/* 00F3C 8099E80C AC400004 */  sw      $zero, 0x0004($v0)         ## 00000004
/* 00F40 8099E810 AC580000 */  sw      $t8, 0x0000($v0)           ## 00000000
/* 00F44 8099E814 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 00F48 8099E818 2409FFFF */  addiu   $t1, $zero, 0xFFFF         ## $t1 = FFFFFFFF
/* 00F4C 8099E81C 3C01809A */  lui     $at, %hi(D_8099EC28)       ## $at = 809A0000
/* 00F50 8099E820 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 00F54 8099E824 AE5902D0 */  sw      $t9, 0x02D0($s2)           ## 000002D0
/* 00F58 8099E828 AC490004 */  sw      $t1, 0x0004($v0)           ## 00000004
/* 00F5C 8099E82C AC480000 */  sw      $t0, 0x0000($v0)           ## 00000000
/* 00F60 8099E830 8FAA00B0 */  lw      $t2, 0x00B0($sp)           
/* 00F64 8099E834 C426EC28 */  lwc1    $f6, %lo(D_8099EC28)($at)  
/* 00F68 8099E838 3C08DB06 */  lui     $t0, 0xDB06                ## $t0 = DB060000
/* 00F6C 8099E83C C5440858 */  lwc1    $f4, 0x0858($t2)           ## 00000858
/* 00F70 8099E840 3C09809A */  lui     $t1, %hi(D_8099EB60)       ## $t1 = 809A0000
/* 00F74 8099E844 26F3024C */  addiu   $s3, $s7, 0x024C           ## $s3 = 0000024C
/* 00F78 8099E848 4604303E */  c.le.s  $f6, $f4                   
/* 00F7C 8099E84C 26F1014C */  addiu   $s1, $s7, 0x014C           ## $s1 = 0000014C
/* 00F80 8099E850 2529EB60 */  addiu   $t1, $t1, %lo(D_8099EB60)  ## $t1 = 8099EB60
/* 00F84 8099E854 35080020 */  ori     $t0, $t0, 0x0020           ## $t0 = DB060020
/* 00F88 8099E858 45000009 */  bc1f    .L8099E880                 
/* 00F8C 8099E85C 0000A025 */  or      $s4, $zero, $zero          ## $s4 = 00000000
/* 00F90 8099E860 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 00F94 8099E864 3C0CFB00 */  lui     $t4, 0xFB00                ## $t4 = FB000000
/* 00F98 8099E868 3C0DFF00 */  lui     $t5, 0xFF00                ## $t5 = FF000000
/* 00F9C 8099E86C 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 00FA0 8099E870 AE4B02D0 */  sw      $t3, 0x02D0($s2)           ## 000002D0
/* 00FA4 8099E874 AC4D0004 */  sw      $t5, 0x0004($v0)           ## 00000004
/* 00FA8 8099E878 10000008 */  beq     $zero, $zero, .L8099E89C   
/* 00FAC 8099E87C AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
.L8099E880:
/* 00FB0 8099E880 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 00FB4 8099E884 3C0FFB00 */  lui     $t7, 0xFB00                ## $t7 = FB000000
/* 00FB8 8099E888 3418FF00 */  ori     $t8, $zero, 0xFF00         ## $t8 = 0000FF00
/* 00FBC 8099E88C 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 00FC0 8099E890 AE4E02D0 */  sw      $t6, 0x02D0($s2)           ## 000002D0
/* 00FC4 8099E894 AC580004 */  sw      $t8, 0x0004($v0)           ## 00000004
/* 00FC8 8099E898 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
.L8099E89C:
/* 00FCC 8099E89C 3C01C1A0 */  lui     $at, 0xC1A0                ## $at = C1A00000
/* 00FD0 8099E8A0 4481E000 */  mtc1    $at, $f28                  ## $f28 = -20.00
/* 00FD4 8099E8A4 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 00FD8 8099E8A8 3C0143A0 */  lui     $at, 0x43A0                ## $at = 43A00000
/* 00FDC 8099E8AC 4481D000 */  mtc1    $at, $f26                  ## $f26 = 320.00
/* 00FE0 8099E8B0 3C01437F */  lui     $at, 0x437F                ## $at = 437F0000
/* 00FE4 8099E8B4 4481C000 */  mtc1    $at, $f24                  ## $f24 = 255.00
/* 00FE8 8099E8B8 24590008 */  addiu   $t9, $v0, 0x0008           ## $t9 = 00000008
/* 00FEC 8099E8BC AE5902D0 */  sw      $t9, 0x02D0($s2)           ## 000002D0
/* 00FF0 8099E8C0 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 00FF4 8099E8C4 4481B000 */  mtc1    $at, $f22                  ## $f22 = 1.00
/* 00FF8 8099E8C8 AC490004 */  sw      $t1, 0x0004($v0)           ## 00000004
/* 00FFC 8099E8CC AC480000 */  sw      $t0, 0x0000($v0)           ## 00000000
.L8099E8D0:
/* 01000 8099E8D0 C6280000 */  lwc1    $f8, 0x0000($s1)           ## 0000014C
/* 01004 8099E8D4 4616403C */  c.lt.s  $f8, $f22                  
/* 01008 8099E8D8 00000000 */  nop
/* 0100C 8099E8DC 45020074 */  bc1fl   .L8099EAB0                 
/* 01010 8099E8E0 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000001
/* 01014 8099E8E4 8E4302D0 */  lw      $v1, 0x02D0($s2)           ## 000002D0
/* 01018 8099E8E8 8FA400B0 */  lw      $a0, 0x00B0($sp)           
/* 0101C 8099E8EC 3C0FFA00 */  lui     $t7, 0xFA00                ## $t7 = FA000000
/* 01020 8099E8F0 246E0008 */  addiu   $t6, $v1, 0x0008           ## $t6 = 00000008
/* 01024 8099E8F4 AE4E02D0 */  sw      $t6, 0x02D0($s2)           ## 000002D0
/* 01028 8099E8F8 AC6F0000 */  sw      $t7, 0x0000($v1)           ## 00000000
/* 0102C 8099E8FC C62A0000 */  lwc1    $f10, 0x0000($s1)          ## 0000014C
/* 01030 8099E900 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 01034 8099E904 3C020403 */  lui     $v0, 0x0403                ## $v0 = 04030000
/* 01038 8099E908 46185402 */  mul.s   $f16, $f10, $f24           
/* 0103C 8099E90C 24427880 */  addiu   $v0, $v0, 0x7880           ## $v0 = 04037880
/* 01040 8099E910 3C0100FF */  lui     $at, 0x00FF                ## $at = 00FF0000
/* 01044 8099E914 3421FFFF */  ori     $at, $at, 0xFFFF           ## $at = 00FFFFFF
/* 01048 8099E918 0041B024 */  and     $s6, $v0, $at              
/* 0104C 8099E91C 00025100 */  sll     $t2, $v0,  4               
/* 01050 8099E920 000A5F02 */  srl     $t3, $t2, 28               
/* 01054 8099E924 4458F800 */  cfc1    $t8, $31
/* 01058 8099E928 44D9F800 */  ctc1    $t9, $31
/* 0105C 8099E92C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 01060 8099E930 3C0D8016 */  lui     $t5, 0x8016                ## $t5 = 80160000
/* 01064 8099E934 460084A4 */  cvt.w.s $f18, $f16                 
/* 01068 8099E938 25AD6FA8 */  addiu   $t5, $t5, 0x6FA8           ## $t5 = 80166FA8
/* 0106C 8099E93C 34211DA0 */  ori     $at, $at, 0x1DA0           ## $at = 00011DA0
/* 01070 8099E940 000B6080 */  sll     $t4, $t3,  2               
/* 01074 8099E944 4459F800 */  cfc1    $t9, $31
/* 01078 8099E948 018DA821 */  addu    $s5, $t4, $t5              
/* 0107C 8099E94C 03C18021 */  addu    $s0, $s8, $at              
/* 01080 8099E950 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 01084 8099E954 13200013 */  beq     $t9, $zero, .L8099E9A4     
/* 01088 8099E958 248409E0 */  addiu   $a0, $a0, 0x09E0           ## $a0 = 000009E0
/* 0108C 8099E95C 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 01090 8099E960 44819000 */  mtc1    $at, $f18                  ## $f18 = 2147483648.00
/* 01094 8099E964 24190001 */  addiu   $t9, $zero, 0x0001         ## $t9 = 00000001
/* 01098 8099E968 46128481 */  sub.s   $f18, $f16, $f18           
/* 0109C 8099E96C 44D9F800 */  ctc1    $t9, $31
/* 010A0 8099E970 00000000 */  nop
/* 010A4 8099E974 460094A4 */  cvt.w.s $f18, $f18                 
/* 010A8 8099E978 4459F800 */  cfc1    $t9, $31
/* 010AC 8099E97C 00000000 */  nop
/* 010B0 8099E980 33390078 */  andi    $t9, $t9, 0x0078           ## $t9 = 00000000
/* 010B4 8099E984 17200005 */  bne     $t9, $zero, .L8099E99C     
/* 010B8 8099E988 00000000 */  nop
/* 010BC 8099E98C 44199000 */  mfc1    $t9, $f18                  
/* 010C0 8099E990 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 010C4 8099E994 10000007 */  beq     $zero, $zero, .L8099E9B4   
/* 010C8 8099E998 0321C825 */  or      $t9, $t9, $at              ## $t9 = 80000000
.L8099E99C:
/* 010CC 8099E99C 10000005 */  beq     $zero, $zero, .L8099E9B4   
/* 010D0 8099E9A0 2419FFFF */  addiu   $t9, $zero, 0xFFFF         ## $t9 = FFFFFFFF
.L8099E9A4:
/* 010D4 8099E9A4 44199000 */  mfc1    $t9, $f18                  
/* 010D8 8099E9A8 00000000 */  nop
/* 010DC 8099E9AC 0720FFFB */  bltz    $t9, .L8099E99C            
/* 010E0 8099E9B0 00000000 */  nop
.L8099E9B4:
/* 010E4 8099E9B4 332800FF */  andi    $t0, $t9, 0x00FF           ## $t0 = 000000FF
/* 010E8 8099E9B8 2401FF00 */  addiu   $at, $zero, 0xFF00         ## $at = FFFFFF00
/* 010EC 8099E9BC 01014825 */  or      $t1, $t0, $at              ## $t1 = FFFFFFFF
/* 010F0 8099E9C0 AC690004 */  sw      $t1, 0x0004($v1)           ## 00000004
/* 010F4 8099E9C4 44D8F800 */  ctc1    $t8, $31
/* 010F8 8099E9C8 C6200000 */  lwc1    $f0, 0x0000($s1)           ## 0000014C
/* 010FC 8099E9CC 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 01100 8099E9D0 46000102 */  mul.s   $f4, $f0, $f0              
/* 01104 8099E9D4 0C03424C */  jal     Matrix_Mult              
/* 01108 8099E9D8 4604B501 */  sub.s   $f20, $f22, $f4            
/* 0110C 8099E9DC C6E20550 */  lwc1    $f2, 0x0550($s7)           ## 00000550
/* 01110 8099E9E0 C6320000 */  lwc1    $f18, 0x0000($s1)          ## 0000014C
/* 01114 8099E9E4 C6700000 */  lwc1    $f16, 0x0000($s3)          ## 0000024C
/* 01118 8099E9E8 46141182 */  mul.s   $f6, $f2, $f20             
/* 0111C 8099E9EC 4602B201 */  sub.s   $f8, $f22, $f2             
/* 01120 8099E9F0 C6640004 */  lwc1    $f4, 0x0004($s3)           ## 00000250
/* 01124 8099E9F4 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 01128 8099E9F8 4612B001 */  sub.s   $f0, $f22, $f18            
/* 0112C 8099E9FC 46083280 */  add.s   $f10, $f6, $f8             
/* 01130 8099EA00 C6680008 */  lwc1    $f8, 0x0008($s3)           ## 00000254
/* 01134 8099EA04 460A8302 */  mul.s   $f12, $f16, $f10           
/* 01138 8099EA08 00000000 */  nop
/* 0113C 8099EA0C 46040182 */  mul.s   $f6, $f0, $f4              
/* 01140 8099EA10 00000000 */  nop
/* 01144 8099EA14 46080402 */  mul.s   $f16, $f0, $f8             
/* 01148 8099EA18 461A3380 */  add.s   $f14, $f6, $f26            
/* 0114C 8099EA1C 461C8280 */  add.s   $f10, $f16, $f28           
/* 01150 8099EA20 44065000 */  mfc1    $a2, $f10                  
/* 01154 8099EA24 0C034261 */  jal     Matrix_Translate              
/* 01158 8099EA28 00000000 */  nop
/* 0115C 8099EA2C C6320000 */  lwc1    $f18, 0x0000($s1)          ## 0000014C
/* 01160 8099EA30 C6E4055C */  lwc1    $f4, 0x055C($s7)           ## 0000055C
/* 01164 8099EA34 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 01168 8099EA38 46049302 */  mul.s   $f12, $f18, $f4            
/* 0116C 8099EA3C 44066000 */  mfc1    $a2, $f12                  
/* 01170 8099EA40 0C0342A3 */  jal     Matrix_Scale              
/* 01174 8099EA44 46006386 */  mov.s   $f14, $f12                 
/* 01178 8099EA48 0C0347F5 */  jal     func_800D1FD4              
/* 0117C 8099EA4C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 01180 8099EA50 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 01184 8099EA54 3C0BDA38 */  lui     $t3, 0xDA38                ## $t3 = DA380000
/* 01188 8099EA58 356B0003 */  ori     $t3, $t3, 0x0003           ## $t3 = DA380003
/* 0118C 8099EA5C 244A0008 */  addiu   $t2, $v0, 0x0008           ## $t2 = 00000008
/* 01190 8099EA60 AE4A02D0 */  sw      $t2, 0x02D0($s2)           ## 000002D0
/* 01194 8099EA64 3C05809A */  lui     $a1, %hi(D_8099EBC0)       ## $a1 = 809A0000
/* 01198 8099EA68 24A5EBC0 */  addiu   $a1, $a1, %lo(D_8099EBC0)  ## $a1 = 8099EBC0
/* 0119C 8099EA6C 02402025 */  or      $a0, $s2, $zero            ## $a0 = 00000000
/* 011A0 8099EA70 240601FA */  addiu   $a2, $zero, 0x01FA         ## $a2 = 000001FA
/* 011A4 8099EA74 AC4B0000 */  sw      $t3, 0x0000($v0)           ## 00000000
/* 011A8 8099EA78 0C0346A2 */  jal     Matrix_NewMtx              
/* 011AC 8099EA7C 00408025 */  or      $s0, $v0, $zero            ## $s0 = 00000000
/* 011B0 8099EA80 AE020004 */  sw      $v0, 0x0004($s0)           ## 00000004
/* 011B4 8099EA84 8E4202D0 */  lw      $v0, 0x02D0($s2)           ## 000002D0
/* 011B8 8099EA88 3C0DDE00 */  lui     $t5, 0xDE00                ## $t5 = DE000000
/* 011BC 8099EA8C 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 011C0 8099EA90 244C0008 */  addiu   $t4, $v0, 0x0008           ## $t4 = 00000008
/* 011C4 8099EA94 AE4C02D0 */  sw      $t4, 0x02D0($s2)           ## 000002D0
/* 011C8 8099EA98 AC4D0000 */  sw      $t5, 0x0000($v0)           ## 00000000
/* 011CC 8099EA9C 8EAE0000 */  lw      $t6, 0x0000($s5)           ## 00000000
/* 011D0 8099EAA0 01D67821 */  addu    $t7, $t6, $s6              
/* 011D4 8099EAA4 01E1C021 */  addu    $t8, $t7, $at              
/* 011D8 8099EAA8 AC580004 */  sw      $t8, 0x0004($v0)           ## 00000004
/* 011DC 8099EAAC 26940001 */  addiu   $s4, $s4, 0x0001           ## $s4 = 00000002
.L8099EAB0:
/* 011E0 8099EAB0 24010040 */  addiu   $at, $zero, 0x0040         ## $at = 00000040
/* 011E4 8099EAB4 2673000C */  addiu   $s3, $s3, 0x000C           ## $s3 = 00000258
/* 011E8 8099EAB8 1681FF85 */  bne     $s4, $at, .L8099E8D0       
/* 011EC 8099EABC 26310004 */  addiu   $s1, $s1, 0x0004           ## $s1 = 00000150
/* 011F0 8099EAC0 3C06809A */  lui     $a2, %hi(D_8099EBD0)       ## $a2 = 809A0000
/* 011F4 8099EAC4 24C6EBD0 */  addiu   $a2, $a2, %lo(D_8099EBD0)  ## $a2 = 8099EBD0
/* 011F8 8099EAC8 27A4009C */  addiu   $a0, $sp, 0x009C           ## $a0 = FFFFFFCC
/* 011FC 8099EACC 02402825 */  or      $a1, $s2, $zero            ## $a1 = 00000000
/* 01200 8099EAD0 0C031AD5 */  jal     Graph_CloseDisp              
/* 01204 8099EAD4 24070203 */  addiu   $a3, $zero, 0x0203         ## $a3 = 00000203
/* 01208 8099EAD8 8FBF0064 */  lw      $ra, 0x0064($sp)           
/* 0120C 8099EADC D7B40018 */  ldc1    $f20, 0x0018($sp)          
/* 01210 8099EAE0 D7B60020 */  ldc1    $f22, 0x0020($sp)          
/* 01214 8099EAE4 D7B80028 */  ldc1    $f24, 0x0028($sp)          
/* 01218 8099EAE8 D7BA0030 */  ldc1    $f26, 0x0030($sp)          
/* 0121C 8099EAEC D7BC0038 */  ldc1    $f28, 0x0038($sp)          
/* 01220 8099EAF0 8FB00040 */  lw      $s0, 0x0040($sp)           
/* 01224 8099EAF4 8FB10044 */  lw      $s1, 0x0044($sp)           
/* 01228 8099EAF8 8FB20048 */  lw      $s2, 0x0048($sp)           
/* 0122C 8099EAFC 8FB3004C */  lw      $s3, 0x004C($sp)           
/* 01230 8099EB00 8FB40050 */  lw      $s4, 0x0050($sp)           
/* 01234 8099EB04 8FB50054 */  lw      $s5, 0x0054($sp)           
/* 01238 8099EB08 8FB60058 */  lw      $s6, 0x0058($sp)           
/* 0123C 8099EB0C 8FB7005C */  lw      $s7, 0x005C($sp)           
/* 01240 8099EB10 8FBE0060 */  lw      $s8, 0x0060($sp)           
/* 01244 8099EB14 03E00008 */  jr      $ra                        
/* 01248 8099EB18 27BD00D0 */  addiu   $sp, $sp, 0x00D0           ## $sp = 00000000


