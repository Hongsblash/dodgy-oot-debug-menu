glabel func_80B388AC
/* 0084C 80B388AC 27BDFF70 */  addiu   $sp, $sp, 0xFF70           ## $sp = FFFFFF70
/* 00850 80B388B0 AFBF008C */  sw      $ra, 0x008C($sp)           
/* 00854 80B388B4 AFBE0088 */  sw      $s8, 0x0088($sp)           
/* 00858 80B388B8 AFB70084 */  sw      $s7, 0x0084($sp)           
/* 0085C 80B388BC AFB60080 */  sw      $s6, 0x0080($sp)           
/* 00860 80B388C0 AFB5007C */  sw      $s5, 0x007C($sp)           
/* 00864 80B388C4 AFB40078 */  sw      $s4, 0x0078($sp)           
/* 00868 80B388C8 AFB30074 */  sw      $s3, 0x0074($sp)           
/* 0086C 80B388CC AFB20070 */  sw      $s2, 0x0070($sp)           
/* 00870 80B388D0 AFB1006C */  sw      $s1, 0x006C($sp)           
/* 00874 80B388D4 AFB00068 */  sw      $s0, 0x0068($sp)           
/* 00878 80B388D8 F7BA0060 */  sdc1    $f26, 0x0060($sp)          
/* 0087C 80B388DC F7B80058 */  sdc1    $f24, 0x0058($sp)          
/* 00880 80B388E0 F7B60050 */  sdc1    $f22, 0x0050($sp)          
/* 00884 80B388E4 F7B40048 */  sdc1    $f20, 0x0048($sp)          
/* 00888 80B388E8 84860158 */  lh      $a2, 0x0158($a0)           ## 00000158
/* 0088C 80B388EC 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00890 80B388F0 00A0F025 */  or      $s8, $a1, $zero            ## $s8 = 00000000
/* 00894 80B388F4 8CB21C44 */  lw      $s2, 0x1C44($a1)           ## 00001C44
/* 00898 80B388F8 8494015E */  lh      $s4, 0x015E($a0)           ## 0000015E
/* 0089C 80B388FC 00008825 */  or      $s1, $zero, $zero          ## $s1 = 00000000
/* 008A0 80B38900 18C00056 */  blez    $a2, .L80B38A5C            
/* 008A4 80B38904 24130001 */  addiu   $s3, $zero, 0x0001         ## $s3 = 00000001
/* 008A8 80B38908 3C014248 */  lui     $at, 0x4248                ## $at = 42480000
/* 008AC 80B3890C 4481D000 */  mtc1    $at, $f26                  ## $f26 = 50.00
/* 008B0 80B38910 3C013F80 */  lui     $at, 0x3F80                ## $at = 3F800000
/* 008B4 80B38914 3C178016 */  lui     $s7, 0x8016                ## $s7 = 80160000
/* 008B8 80B38918 3C1580B4 */  lui     $s5, %hi(D_80B39080)       ## $s5 = 80B40000
/* 008BC 80B3891C 4481C000 */  mtc1    $at, $f24                  ## $f24 = 1.00
/* 008C0 80B38920 26B59080 */  addiu   $s5, $s5, %lo(D_80B39080)  ## $s5 = 80B39080
/* 008C4 80B38924 26F7FA90 */  addiu   $s7, $s7, 0xFA90           ## $s7 = 8015FA90
/* 008C8 80B38928 2416000C */  addiu   $s6, $zero, 0x000C         ## $s6 = 0000000C
.L80B3892C:
/* 008CC 80B3892C 02931824 */  and     $v1, $s4, $s3              
/* 008D0 80B38930 54600047 */  bnel    $v1, $zero, .L80B38A50     
/* 008D4 80B38934 26310001 */  addiu   $s1, $s1, 0x0001           ## $s1 = 00000001
/* 008D8 80B38938 02360019 */  multu   $s1, $s6                   
/* 008DC 80B3893C C6440024 */  lwc1    $f4, 0x0024($s2)           ## 00000024
/* 008E0 80B38940 C6460028 */  lwc1    $f6, 0x0028($s2)           ## 00000028
/* 008E4 80B38944 C648002C */  lwc1    $f8, 0x002C($s2)           ## 0000002C
/* 008E8 80B38948 241800FF */  addiu   $t8, $zero, 0x00FF         ## $t8 = 000000FF
/* 008EC 80B3894C 241900FF */  addiu   $t9, $zero, 0x00FF         ## $t9 = 000000FF
/* 008F0 80B38950 24080004 */  addiu   $t0, $zero, 0x0004         ## $t0 = 00000004
/* 008F4 80B38954 00007012 */  mflo    $t6                        
/* 008F8 80B38958 02AE1021 */  addu    $v0, $s5, $t6              
/* 008FC 80B3895C C4540000 */  lwc1    $f20, 0x0000($v0)          ## 00000000
/* 00900 80B38960 C4560004 */  lwc1    $f22, 0x0004($v0)          ## 00000004
/* 00904 80B38964 C4500008 */  lwc1    $f16, 0x0008($v0)          ## 00000008
/* 00908 80B38968 46142081 */  sub.s   $f2, $f4, $f20             
/* 0090C 80B3896C 46163301 */  sub.s   $f12, $f6, $f22            
/* 00910 80B38970 46021282 */  mul.s   $f10, $f2, $f2             
/* 00914 80B38974 46104381 */  sub.s   $f14, $f8, $f16            
/* 00918 80B38978 460C6482 */  mul.s   $f18, $f12, $f12           
/* 0091C 80B3897C 46125100 */  add.s   $f4, $f10, $f18            
/* 00920 80B38980 460E7182 */  mul.s   $f6, $f14, $f14            
/* 00924 80B38984 46062000 */  add.s   $f0, $f4, $f6              
/* 00928 80B38988 46000004 */  sqrt.s  $f0, $f0                   
/* 0092C 80B3898C 461A003C */  c.lt.s  $f0, $f26                  
/* 00930 80B38990 00000000 */  nop
/* 00934 80B38994 45020016 */  bc1fl   .L80B389F0                 
/* 00938 80B38998 8EEC0000 */  lw      $t4, 0x0000($s7)           ## 8015FA90
/* 0093C 80B3899C 5460003F */  bnel    $v1, $zero, .L80B38A9C     
/* 00940 80B389A0 8FBF008C */  lw      $ra, 0x008C($sp)           
/* 00944 80B389A4 86020168 */  lh      $v0, 0x0168($s0)           ## 00000168
/* 00948 80B389A8 1622000C */  bne     $s1, $v0, .L80B389DC       
/* 0094C 80B389AC 24490001 */  addiu   $t1, $v0, 0x0001           ## $t1 = 00000001
/* 00950 80B389B0 860F015E */  lh      $t7, 0x015E($s0)           ## 0000015E
/* 00954 80B389B4 86190160 */  lh      $t9, 0x0160($s0)           ## 00000160
/* 00958 80B389B8 860A016A */  lh      $t2, 0x016A($s0)           ## 0000016A
/* 0095C 80B389BC 01F3C025 */  or      $t8, $t7, $s3              ## $t8 = 00000001
/* 00960 80B389C0 27280001 */  addiu   $t0, $t9, 0x0001           ## $t0 = 00000100
/* 00964 80B389C4 254B0051 */  addiu   $t3, $t2, 0x0051           ## $t3 = 00000051
/* 00968 80B389C8 A618015E */  sh      $t8, 0x015E($s0)           ## 0000015E
/* 0096C 80B389CC A6080160 */  sh      $t0, 0x0160($s0)           ## 00000160
/* 00970 80B389D0 A6090168 */  sh      $t1, 0x0168($s0)           ## 00000168
/* 00974 80B389D4 10000030 */  beq     $zero, $zero, .L80B38A98   
/* 00978 80B389D8 A60B015C */  sh      $t3, 0x015C($s0)           ## 0000015C
.L80B389DC:
/* 0097C 80B389DC 0C00B55C */  jal     Actor_Kill
              
/* 00980 80B389E0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00984 80B389E4 1000002D */  beq     $zero, $zero, .L80B38A9C   
/* 00988 80B389E8 8FBF008C */  lw      $ra, 0x008C($sp)           
/* 0098C 80B389EC 8EEC0000 */  lw      $t4, 0x0000($s7)           ## 8015FA90
.L80B389F0:
/* 00990 80B389F0 858D12D4 */  lh      $t5, 0x12D4($t4)           ## 000012D4
/* 00994 80B389F4 51A00016 */  beql    $t5, $zero, .L80B38A50     
/* 00998 80B389F8 26310001 */  addiu   $s1, $s1, 0x0001           ## $s1 = 00000002
/* 0099C 80B389FC 860E0032 */  lh      $t6, 0x0032($s0)           ## 00000032
/* 009A0 80B38A00 86070030 */  lh      $a3, 0x0030($s0)           ## 00000030
/* 009A4 80B38A04 44068000 */  mfc1    $a2, $f16                  
/* 009A8 80B38A08 AFAE0010 */  sw      $t6, 0x0010($sp)           
/* 009AC 80B38A0C 860F0034 */  lh      $t7, 0x0034($s0)           ## 00000034
/* 009B0 80B38A10 AFA80034 */  sw      $t0, 0x0034($sp)           
/* 009B4 80B38A14 AFB90030 */  sw      $t9, 0x0030($sp)           
/* 009B8 80B38A18 AFB8002C */  sw      $t8, 0x002C($sp)           
/* 009BC 80B38A1C AFA00028 */  sw      $zero, 0x0028($sp)         
/* 009C0 80B38A20 AFA00024 */  sw      $zero, 0x0024($sp)         
/* 009C4 80B38A24 E7B80020 */  swc1    $f24, 0x0020($sp)          
/* 009C8 80B38A28 E7B8001C */  swc1    $f24, 0x001C($sp)          
/* 009CC 80B38A2C E7B80018 */  swc1    $f24, 0x0018($sp)          
/* 009D0 80B38A30 AFAF0014 */  sw      $t7, 0x0014($sp)           
/* 009D4 80B38A34 8FC90000 */  lw      $t1, 0x0000($s8)           ## 00000000
/* 009D8 80B38A38 4600A306 */  mov.s   $f12, $f20                 
/* 009DC 80B38A3C 4600B386 */  mov.s   $f14, $f22                 
/* 009E0 80B38A40 0C018FA7 */  jal     DebugDisplay_AddObject
              
/* 009E4 80B38A44 AFA90038 */  sw      $t1, 0x0038($sp)           
/* 009E8 80B38A48 86060158 */  lh      $a2, 0x0158($s0)           ## 00000158
/* 009EC 80B38A4C 26310001 */  addiu   $s1, $s1, 0x0001           ## $s1 = 00000003
.L80B38A50:
/* 009F0 80B38A50 0226082A */  slt     $at, $s1, $a2              
/* 009F4 80B38A54 1420FFB5 */  bne     $at, $zero, .L80B3892C     
/* 009F8 80B38A58 00139840 */  sll     $s3, $s3,  1               
.L80B38A5C:
/* 009FC 80B38A5C 860A015C */  lh      $t2, 0x015C($s0)           ## 0000015C
/* 00A00 80B38A60 24010001 */  addiu   $at, $zero, 0x0001         ## $at = 00000001
/* 00A04 80B38A64 55410006 */  bnel    $t2, $at, .L80B38A80       
/* 00A08 80B38A68 860B0160 */  lh      $t3, 0x0160($s0)           ## 00000160
/* 00A0C 80B38A6C 0C00B55C */  jal     Actor_Kill
              
/* 00A10 80B38A70 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00A14 80B38A74 10000009 */  beq     $zero, $zero, .L80B38A9C   
/* 00A18 80B38A78 8FBF008C */  lw      $ra, 0x008C($sp)           
/* 00A1C 80B38A7C 860B0160 */  lh      $t3, 0x0160($s0)           ## 00000160
.L80B38A80:
/* 00A20 80B38A80 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00A24 80B38A84 03C02825 */  or      $a1, $s8, $zero            ## $a1 = 00000000
/* 00A28 80B38A88 54CB0004 */  bnel    $a2, $t3, .L80B38A9C       
/* 00A2C 80B38A8C 8FBF008C */  lw      $ra, 0x008C($sp)           
/* 00A30 80B38A90 0C2CE028 */  jal     func_80B380A0              
/* 00A34 80B38A94 24060001 */  addiu   $a2, $zero, 0x0001         ## $a2 = 00000001
.L80B38A98:
/* 00A38 80B38A98 8FBF008C */  lw      $ra, 0x008C($sp)           
.L80B38A9C:
/* 00A3C 80B38A9C D7B40048 */  ldc1    $f20, 0x0048($sp)          
/* 00A40 80B38AA0 D7B60050 */  ldc1    $f22, 0x0050($sp)          
/* 00A44 80B38AA4 D7B80058 */  ldc1    $f24, 0x0058($sp)          
/* 00A48 80B38AA8 D7BA0060 */  ldc1    $f26, 0x0060($sp)          
/* 00A4C 80B38AAC 8FB00068 */  lw      $s0, 0x0068($sp)           
/* 00A50 80B38AB0 8FB1006C */  lw      $s1, 0x006C($sp)           
/* 00A54 80B38AB4 8FB20070 */  lw      $s2, 0x0070($sp)           
/* 00A58 80B38AB8 8FB30074 */  lw      $s3, 0x0074($sp)           
/* 00A5C 80B38ABC 8FB40078 */  lw      $s4, 0x0078($sp)           
/* 00A60 80B38AC0 8FB5007C */  lw      $s5, 0x007C($sp)           
/* 00A64 80B38AC4 8FB60080 */  lw      $s6, 0x0080($sp)           
/* 00A68 80B38AC8 8FB70084 */  lw      $s7, 0x0084($sp)           
/* 00A6C 80B38ACC 8FBE0088 */  lw      $s8, 0x0088($sp)           
/* 00A70 80B38AD0 03E00008 */  jr      $ra                        
/* 00A74 80B38AD4 27BD0090 */  addiu   $sp, $sp, 0x0090           ## $sp = 00000000


