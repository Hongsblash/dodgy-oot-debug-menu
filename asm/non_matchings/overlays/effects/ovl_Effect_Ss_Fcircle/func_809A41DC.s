glabel func_809A41DC
/* 000BC 809A41DC 27BDFF80 */  addiu   $sp, $sp, 0xFF80           ## $sp = FFFFFF80
/* 000C0 809A41E0 AFB10038 */  sw      $s1, 0x0038($sp)           
/* 000C4 809A41E4 AFBF003C */  sw      $ra, 0x003C($sp)           
/* 000C8 809A41E8 AFB00034 */  sw      $s0, 0x0034($sp)           
/* 000CC 809A41EC AFA40080 */  sw      $a0, 0x0080($sp)           
/* 000D0 809A41F0 AFA50084 */  sw      $a1, 0x0084($sp)           
/* 000D4 809A41F4 8C900000 */  lw      $s0, 0x0000($a0)           ## 00000000
/* 000D8 809A41F8 00C08825 */  or      $s1, $a2, $zero            ## $s1 = 00000000
/* 000DC 809A41FC 3C06809A */  lui     $a2, %hi(D_809A45A0)       ## $a2 = 809A0000
/* 000E0 809A4200 24C645A0 */  addiu   $a2, $a2, %lo(D_809A45A0)  ## $a2 = 809A45A0
/* 000E4 809A4204 27A40058 */  addiu   $a0, $sp, 0x0058           ## $a0 = FFFFFFD8
/* 000E8 809A4208 24070095 */  addiu   $a3, $zero, 0x0095         ## $a3 = 00000095
/* 000EC 809A420C 0C031AB1 */  jal     func_800C6AC4              
/* 000F0 809A4210 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 000F4 809A4214 862F005C */  lh      $t7, 0x005C($s1)           ## 0000005C
/* 000F8 809A4218 3C01809A */  lui     $at, %hi(D_809A45DC)       ## $at = 809A0000
/* 000FC 809A421C C42245DC */  lwc1    $f2, %lo(D_809A45DC)($at)  
/* 00100 809A4220 448F2000 */  mtc1    $t7, $f4                   ## $f4 = 0.00
/* 00104 809A4224 3C01809A */  lui     $at, %hi(D_809A45E0)       ## $at = 809A0000
/* 00108 809A4228 C42845E0 */  lwc1    $f8, %lo(D_809A45E0)($at)  
/* 0010C 809A422C 468021A0 */  cvt.s.w $f6, $f4                   
/* 00110 809A4230 86380056 */  lh      $t8, 0x0056($s1)           ## 00000056
/* 00114 809A4234 3C013F00 */  lui     $at, 0x3F00                ## $at = 3F000000
/* 00118 809A4238 44818000 */  mtc1    $at, $f16                  ## $f16 = 0.50
/* 0011C 809A423C 44982000 */  mtc1    $t8, $f4                   ## $f4 = 0.00
/* 00120 809A4240 86390052 */  lh      $t9, 0x0052($s1)           ## 00000052
/* 00124 809A4244 46083282 */  mul.s   $f10, $f6, $f8             
/* 00128 809A4248 3C01809A */  lui     $at, %hi(D_809A45E4)       ## $at = 809A0000
/* 0012C 809A424C 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 00130 809A4250 468021A0 */  cvt.s.w $f6, $f4                   
/* 00134 809A4254 460A8480 */  add.s   $f18, $f16, $f10           
/* 00138 809A4258 44995000 */  mtc1    $t9, $f10                  ## $f10 = 0.00
/* 0013C 809A425C C43045E4 */  lwc1    $f16, %lo(D_809A45E4)($at) 
/* 00140 809A4260 46123202 */  mul.s   $f8, $f6, $f18             
/* 00144 809A4264 46805120 */  cvt.s.w $f4, $f10                  
/* 00148 809A4268 46104002 */  mul.s   $f0, $f8, $f16             
/* 0014C 809A426C 00000000 */  nop
/* 00150 809A4270 46022182 */  mul.s   $f6, $f4, $f2              
/* 00154 809A4274 00000000 */  nop
/* 00158 809A4278 46003482 */  mul.s   $f18, $f6, $f0             
/* 0015C 809A427C E7B20074 */  swc1    $f18, 0x0074($sp)          
/* 00160 809A4280 86290050 */  lh      $t1, 0x0050($s1)           ## 00000050
/* 00164 809A4284 44894000 */  mtc1    $t1, $f8                   ## $f8 = 0.00
/* 00168 809A4288 00000000 */  nop
/* 0016C 809A428C 46804420 */  cvt.s.w $f16, $f8                  
/* 00170 809A4290 46028282 */  mul.s   $f10, $f16, $f2            
/* 00174 809A4294 00000000 */  nop
/* 00178 809A4298 46005102 */  mul.s   $f4, $f10, $f0             
/* 0017C 809A429C E7A40070 */  swc1    $f4, 0x0070($sp)           
/* 00180 809A42A0 8E260008 */  lw      $a2, 0x0008($s1)           ## 00000008
/* 00184 809A42A4 C62E0004 */  lwc1    $f14, 0x0004($s1)          ## 00000004
/* 00188 809A42A8 0C034261 */  jal     Matrix_Translate              
/* 0018C 809A42AC C62C0000 */  lwc1    $f12, 0x0000($s1)          ## 00000000
/* 00190 809A42B0 C7AC0070 */  lwc1    $f12, 0x0070($sp)          
/* 00194 809A42B4 C7AE0074 */  lwc1    $f14, 0x0074($sp)          
/* 00198 809A42B8 24070001 */  addiu   $a3, $zero, 0x0001         ## $a3 = 00000001
/* 0019C 809A42BC 44066000 */  mfc1    $a2, $f12                  
/* 001A0 809A42C0 0C0342A3 */  jal     Matrix_Scale              
/* 001A4 809A42C4 00000000 */  nop
/* 001A8 809A42C8 862A0054 */  lh      $t2, 0x0054($s1)           ## 00000054
/* 001AC 809A42CC 3C01809A */  lui     $at, %hi(D_809A45E8)       ## $at = 809A0000
/* 001B0 809A42D0 C42845E8 */  lwc1    $f8, %lo(D_809A45E8)($at)  
/* 001B4 809A42D4 448A3000 */  mtc1    $t2, $f6                   ## $f6 = 0.00
/* 001B8 809A42D8 24050001 */  addiu   $a1, $zero, 0x0001         ## $a1 = 00000001
/* 001BC 809A42DC 468034A0 */  cvt.s.w $f18, $f6                  
/* 001C0 809A42E0 46089302 */  mul.s   $f12, $f18, $f8            
/* 001C4 809A42E4 0C034348 */  jal     Matrix_RotateY              
/* 001C8 809A42E8 00000000 */  nop
/* 001CC 809A42EC 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 001D0 809A42F0 3C0CDA38 */  lui     $t4, 0xDA38                ## $t4 = DA380000
/* 001D4 809A42F4 358C0003 */  ori     $t4, $t4, 0x0003           ## $t4 = DA380003
/* 001D8 809A42F8 244B0008 */  addiu   $t3, $v0, 0x0008           ## $t3 = 00000008
/* 001DC 809A42FC AE0B02D0 */  sw      $t3, 0x02D0($s0)           ## 000002D0
/* 001E0 809A4300 3C05809A */  lui     $a1, %hi(D_809A45B4)       ## $a1 = 809A0000
/* 001E4 809A4304 AC4C0000 */  sw      $t4, 0x0000($v0)           ## 00000000
/* 001E8 809A4308 24A545B4 */  addiu   $a1, $a1, %lo(D_809A45B4)  ## $a1 = 809A45B4
/* 001EC 809A430C 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 001F0 809A4310 240600A3 */  addiu   $a2, $zero, 0x00A3         ## $a2 = 000000A3
/* 001F4 809A4314 0C0346A2 */  jal     Matrix_NewMtx              
/* 001F8 809A4318 AFA20054 */  sw      $v0, 0x0054($sp)           
/* 001FC 809A431C 8FA30054 */  lw      $v1, 0x0054($sp)           
/* 00200 809A4320 AC620004 */  sw      $v0, 0x0004($v1)           ## 00000004
/* 00204 809A4324 8FAD0080 */  lw      $t5, 0x0080($sp)           
/* 00208 809A4328 0C024F61 */  jal     func_80093D84              
/* 0020C 809A432C 8DA40000 */  lw      $a0, 0x0000($t5)           ## 00000000
/* 00210 809A4330 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00214 809A4334 3C0FDB06 */  lui     $t7, 0xDB06                ## $t7 = DB060000
/* 00218 809A4338 35EF0020 */  ori     $t7, $t7, 0x0020           ## $t7 = DB060020
/* 0021C 809A433C 244E0008 */  addiu   $t6, $v0, 0x0008           ## $t6 = 00000008
/* 00220 809A4340 AE0E02D0 */  sw      $t6, 0x02D0($s0)           ## 000002D0
/* 00224 809A4344 AC4F0000 */  sw      $t7, 0x0000($v0)           ## 00000000
/* 00228 809A4348 8FB80080 */  lw      $t8, 0x0080($sp)           
/* 0022C 809A434C 3C030001 */  lui     $v1, 0x0001                ## $v1 = 00010000
/* 00230 809A4350 240E0040 */  addiu   $t6, $zero, 0x0040         ## $t6 = 00000040
/* 00234 809A4354 00781821 */  addu    $v1, $v1, $t8              
/* 00238 809A4358 8C631DE4 */  lw      $v1, 0x1DE4($v1)           ## 00011DE4
/* 0023C 809A435C 8F040000 */  lw      $a0, 0x0000($t8)           ## 00000000
/* 00240 809A4360 24190020 */  addiu   $t9, $zero, 0x0020         ## $t9 = 00000020
/* 00244 809A4364 00030823 */  subu    $at, $zero, $v1            
/* 00248 809A4368 00015900 */  sll     $t3, $at,  4               
/* 0024C 809A436C 01615823 */  subu    $t3, $t3, $at              
/* 00250 809A4370 316C00FF */  andi    $t4, $t3, 0x00FF           ## $t4 = 00000000
/* 00254 809A4374 24090040 */  addiu   $t1, $zero, 0x0040         ## $t1 = 00000040
/* 00258 809A4378 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 0025C 809A437C 240D0020 */  addiu   $t5, $zero, 0x0020         ## $t5 = 00000020
/* 00260 809A4380 AFAD0024 */  sw      $t5, 0x0024($sp)           
/* 00264 809A4384 AFAA0018 */  sw      $t2, 0x0018($sp)           
/* 00268 809A4388 AFA90014 */  sw      $t1, 0x0014($sp)           
/* 0026C 809A438C AFAC0020 */  sw      $t4, 0x0020($sp)           
/* 00270 809A4390 AFB90010 */  sw      $t9, 0x0010($sp)           
/* 00274 809A4394 AFAE0028 */  sw      $t6, 0x0028($sp)           
/* 00278 809A4398 AFA0001C */  sw      $zero, 0x001C($sp)         
/* 0027C 809A439C 00002825 */  or      $a1, $zero, $zero          ## $a1 = 00000000
/* 00280 809A43A0 00003825 */  or      $a3, $zero, $zero          ## $a3 = 00000000
/* 00284 809A43A4 AFA20050 */  sw      $v0, 0x0050($sp)           
/* 00288 809A43A8 0C0253D0 */  jal     Draw_TwoTexScroll              
/* 0028C 809A43AC 3066007F */  andi    $a2, $v1, 0x007F           ## $a2 = 00000000
/* 00290 809A43B0 8FA80050 */  lw      $t0, 0x0050($sp)           
/* 00294 809A43B4 3C18FA00 */  lui     $t8, 0xFA00                ## $t8 = FA000000
/* 00298 809A43B8 37188080 */  ori     $t8, $t8, 0x8080           ## $t8 = FA008080
/* 0029C 809A43BC AD020004 */  sw      $v0, 0x0004($t0)           ## 00000004
/* 002A0 809A43C0 8E0302D0 */  lw      $v1, 0x02D0($s0)           ## 000002D0
/* 002A4 809A43C4 3C01414C */  lui     $at, 0x414C                ## $at = 414C0000
/* 002A8 809A43C8 44812000 */  mtc1    $at, $f4                   ## $f4 = 12.75
/* 002AC 809A43CC 246F0008 */  addiu   $t7, $v1, 0x0008           ## $t7 = 00000008
/* 002B0 809A43D0 AE0F02D0 */  sw      $t7, 0x02D0($s0)           ## 000002D0
/* 002B4 809A43D4 AC780000 */  sw      $t8, 0x0000($v1)           ## 00000000
/* 002B8 809A43D8 8639005C */  lh      $t9, 0x005C($s1)           ## 0000005C
/* 002BC 809A43DC 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 002C0 809A43E0 3C014F00 */  lui     $at, 0x4F00                ## $at = 4F000000
/* 002C4 809A43E4 44998000 */  mtc1    $t9, $f16                  ## $f16 = 0.00
/* 002C8 809A43E8 3C0FFF00 */  lui     $t7, 0xFF00                ## $t7 = FF000000
/* 002CC 809A43EC 3C0EFB00 */  lui     $t6, 0xFB00                ## $t6 = FB000000
/* 002D0 809A43F0 468082A0 */  cvt.s.w $f10, $f16                 
/* 002D4 809A43F4 3C19DE00 */  lui     $t9, 0xDE00                ## $t9 = DE000000
/* 002D8 809A43F8 27A40058 */  addiu   $a0, $sp, 0x0058           ## $a0 = FFFFFFD8
/* 002DC 809A43FC 46045182 */  mul.s   $f6, $f10, $f4             
/* 002E0 809A4400 4449F800 */  cfc1    $t1, $31
/* 002E4 809A4404 44CAF800 */  ctc1    $t2, $31
/* 002E8 809A4408 00000000 */  nop
/* 002EC 809A440C 460034A4 */  cvt.w.s $f18, $f6                  
/* 002F0 809A4410 444AF800 */  cfc1    $t2, $31
/* 002F4 809A4414 00000000 */  nop
/* 002F8 809A4418 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 002FC 809A441C 51400013 */  beql    $t2, $zero, .L809A446C     
/* 00300 809A4420 440A9000 */  mfc1    $t2, $f18                  
/* 00304 809A4424 44819000 */  mtc1    $at, $f18                  ## $f18 = 2147483648.00
/* 00308 809A4428 240A0001 */  addiu   $t2, $zero, 0x0001         ## $t2 = 00000001
/* 0030C 809A442C 46123481 */  sub.s   $f18, $f6, $f18            
/* 00310 809A4430 44CAF800 */  ctc1    $t2, $31
/* 00314 809A4434 00000000 */  nop
/* 00318 809A4438 460094A4 */  cvt.w.s $f18, $f18                 
/* 0031C 809A443C 444AF800 */  cfc1    $t2, $31
/* 00320 809A4440 00000000 */  nop
/* 00324 809A4444 314A0078 */  andi    $t2, $t2, 0x0078           ## $t2 = 00000000
/* 00328 809A4448 15400005 */  bne     $t2, $zero, .L809A4460     
/* 0032C 809A444C 00000000 */  nop
/* 00330 809A4450 440A9000 */  mfc1    $t2, $f18                  
/* 00334 809A4454 3C018000 */  lui     $at, 0x8000                ## $at = 80000000
/* 00338 809A4458 10000007 */  beq     $zero, $zero, .L809A4478   
/* 0033C 809A445C 01415025 */  or      $t2, $t2, $at              ## $t2 = 80000000
.L809A4460:
/* 00340 809A4460 10000005 */  beq     $zero, $zero, .L809A4478   
/* 00344 809A4464 240AFFFF */  addiu   $t2, $zero, 0xFFFF         ## $t2 = FFFFFFFF
/* 00348 809A4468 440A9000 */  mfc1    $t2, $f18                  
.L809A446C:
/* 0034C 809A446C 00000000 */  nop
/* 00350 809A4470 0540FFFB */  bltz    $t2, .L809A4460            
/* 00354 809A4474 00000000 */  nop
.L809A4478:
/* 00358 809A4478 314B00FF */  andi    $t3, $t2, 0x00FF           ## $t3 = 000000FF
/* 0035C 809A447C 3C01FFDC */  lui     $at, 0xFFDC                ## $at = FFDC0000
/* 00360 809A4480 01616025 */  or      $t4, $t3, $at              ## $t4 = FFDC00FF
/* 00364 809A4484 AC6C0004 */  sw      $t4, 0x0004($v1)           ## 00000004
/* 00368 809A4488 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 0036C 809A448C 44C9F800 */  ctc1    $t1, $31
/* 00370 809A4490 3C06809A */  lui     $a2, %hi(D_809A45C8)       ## $a2 = 809A0000
/* 00374 809A4494 244D0008 */  addiu   $t5, $v0, 0x0008           ## $t5 = 00000008
/* 00378 809A4498 AE0D02D0 */  sw      $t5, 0x02D0($s0)           ## 000002D0
/* 0037C 809A449C AC4F0004 */  sw      $t7, 0x0004($v0)           ## 00000004
/* 00380 809A44A0 AC4E0000 */  sw      $t6, 0x0000($v0)           ## 00000000
/* 00384 809A44A4 8E0202D0 */  lw      $v0, 0x02D0($s0)           ## 000002D0
/* 00388 809A44A8 24C645C8 */  addiu   $a2, $a2, %lo(D_809A45C8)  ## $a2 = 809A45C8
/* 0038C 809A44AC 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 00390 809A44B0 24580008 */  addiu   $t8, $v0, 0x0008           ## $t8 = 00000008
/* 00394 809A44B4 AE1802D0 */  sw      $t8, 0x02D0($s0)           ## 000002D0
/* 00398 809A44B8 AC590000 */  sw      $t9, 0x0000($v0)           ## 00000000
/* 0039C 809A44BC 8E290038 */  lw      $t1, 0x0038($s1)           ## 00000038
/* 003A0 809A44C0 240700BA */  addiu   $a3, $zero, 0x00BA         ## $a3 = 000000BA
/* 003A4 809A44C4 0C031AD5 */  jal     func_800C6B54              
/* 003A8 809A44C8 AC490004 */  sw      $t1, 0x0004($v0)           ## 00000004
/* 003AC 809A44CC 8FBF003C */  lw      $ra, 0x003C($sp)           
/* 003B0 809A44D0 8FB00034 */  lw      $s0, 0x0034($sp)           
/* 003B4 809A44D4 8FB10038 */  lw      $s1, 0x0038($sp)           
/* 003B8 809A44D8 03E00008 */  jr      $ra                        
/* 003BC 809A44DC 27BD0080 */  addiu   $sp, $sp, 0x0080           ## $sp = 00000000


