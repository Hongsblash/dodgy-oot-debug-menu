glabel func_80897B48
/* 001D8 80897B48 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 001DC 80897B4C 3C018089 */  lui     $at, %hi(D_80897FC4)       ## $at = 80890000
/* 001E0 80897B50 C4207FC4 */  lwc1    $f0, %lo(D_80897FC4)($at)  
/* 001E4 80897B54 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 001E8 80897B58 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 001EC 80897B5C AFA5002C */  sw      $a1, 0x002C($sp)           
/* 001F0 80897B60 3C01C2C8 */  lui     $at, 0xC2C8                ## $at = C2C80000
/* 001F4 80897B64 44812000 */  mtc1    $at, $f4                   ## $f4 = -100.00
/* 001F8 80897B68 C4860028 */  lwc1    $f6, 0x0028($a0)           ## 00000028
/* 001FC 80897B6C 3C014020 */  lui     $at, 0x4020                ## $at = 40200000
/* 00200 80897B70 44815000 */  mtc1    $at, $f10                  ## $f10 = 2.50
/* 00204 80897B74 46062201 */  sub.s   $f8, $f4, $f6              
/* 00208 80897B78 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 0020C 80897B7C 460A4082 */  mul.s   $f2, $f8, $f10             
/* 00210 80897B80 4600103C */  c.lt.s  $f2, $f0                   
/* 00214 80897B84 00000000 */  nop
/* 00218 80897B88 45020003 */  bc1fl   .L80897B98                 
/* 0021C 80897B8C 46001004 */  sqrt.s  $f0, $f2                   
/* 00220 80897B90 46000086 */  mov.s   $f2, $f0                   
/* 00224 80897B94 46001004 */  sqrt.s  $f0, $f2                   
.L80897B98:
/* 00228 80897B98 C61001B0 */  lwc1    $f16, 0x01B0($s0)          ## 000001B0
/* 0022C 80897B9C 86040032 */  lh      $a0, 0x0032($s0)           ## 00000032
/* 00230 80897BA0 46100482 */  mul.s   $f18, $f0, $f16            
/* 00234 80897BA4 0C01DE1C */  jal     Math_Sins
              ## sins?
/* 00238 80897BA8 E6120068 */  swc1    $f18, 0x0068($s0)          ## 00000068
/* 0023C 80897BAC C6040068 */  lwc1    $f4, 0x0068($s0)           ## 00000068
/* 00240 80897BB0 86040032 */  lh      $a0, 0x0032($s0)           ## 00000032
/* 00244 80897BB4 46040182 */  mul.s   $f6, $f0, $f4              
/* 00248 80897BB8 0C01DE0D */  jal     Math_Coss
              ## coss?
/* 0024C 80897BBC E606005C */  swc1    $f6, 0x005C($s0)           ## 0000005C
/* 00250 80897BC0 C6080068 */  lwc1    $f8, 0x0068($s0)           ## 00000068
/* 00254 80897BC4 C60A0024 */  lwc1    $f10, 0x0024($s0)          ## 00000024
/* 00258 80897BC8 C610005C */  lwc1    $f16, 0x005C($s0)          ## 0000005C
/* 0025C 80897BCC 46080302 */  mul.s   $f12, $f0, $f8             
/* 00260 80897BD0 C604002C */  lwc1    $f4, 0x002C($s0)           ## 0000002C
/* 00264 80897BD4 46105480 */  add.s   $f18, $f10, $f16           
/* 00268 80897BD8 3C018089 */  lui     $at, %hi(D_80897FC8)       ## $at = 80890000
/* 0026C 80897BDC 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 00270 80897BE0 3C064000 */  lui     $a2, 0x4000                ## $a2 = 40000000
/* 00274 80897BE4 E6120024 */  swc1    $f18, 0x0024($s0)          ## 00000024
/* 00278 80897BE8 460C2180 */  add.s   $f6, $f4, $f12             
/* 0027C 80897BEC E60C0064 */  swc1    $f12, 0x0064($s0)          ## 00000064
/* 00280 80897BF0 C6020024 */  lwc1    $f2, 0x0024($s0)           ## 00000024
/* 00284 80897BF4 E606002C */  swc1    $f6, 0x002C($s0)           ## 0000002C
/* 00288 80897BF8 C4287FC8 */  lwc1    $f8, %lo(D_80897FC8)($at)  
/* 0028C 80897BFC 3C018089 */  lui     $at, %hi(D_80897FCC)       ## $at = 80890000
/* 00290 80897C00 4602403C */  c.lt.s  $f8, $f2                   
/* 00294 80897C04 00000000 */  nop
/* 00298 80897C08 4500000A */  bc1f    .L80897C34                 
/* 0029C 80897C0C 00000000 */  nop
/* 002A0 80897C10 C42A7FCC */  lwc1    $f10, %lo(D_80897FCC)($at) 
/* 002A4 80897C14 3C018089 */  lui     $at, %hi(D_80897FD0)       ## $at = 80890000
/* 002A8 80897C18 460A103C */  c.lt.s  $f2, $f10                  
/* 002AC 80897C1C 00000000 */  nop
/* 002B0 80897C20 45000004 */  bc1f    .L80897C34                 
/* 002B4 80897C24 00000000 */  nop
/* 002B8 80897C28 C4307FD0 */  lwc1    $f16, %lo(D_80897FD0)($at) 
/* 002BC 80897C2C 1000000F */  beq     $zero, $zero, .L80897C6C   
/* 002C0 80897C30 E6100028 */  swc1    $f16, 0x0028($s0)          ## 00000028
.L80897C34:
/* 002C4 80897C34 3C018089 */  lui     $at, %hi(D_80897FD4)       ## $at = 80890000
/* 002C8 80897C38 C4327FD4 */  lwc1    $f18, %lo(D_80897FD4)($at) 
/* 002CC 80897C3C 3C018089 */  lui     $at, %hi(D_80897FD8)       ## $at = 80890000
/* 002D0 80897C40 C4247FD8 */  lwc1    $f4, %lo(D_80897FD8)($at)  
/* 002D4 80897C44 46029301 */  sub.s   $f12, $f18, $f2            
/* 002D8 80897C48 3C0142CE */  lui     $at, 0x42CE                ## $at = 42CE0000
/* 002DC 80897C4C 44813000 */  mtc1    $at, $f6                   ## $f6 = 103.00
/* 002E0 80897C50 3C018089 */  lui     $at, %hi(D_80897FDC)       ## $at = 80890000
/* 002E4 80897C54 46006005 */  abs.s   $f0, $f12                  
/* 002E8 80897C58 C4307FDC */  lwc1    $f16, %lo(D_80897FDC)($at) 
/* 002EC 80897C5C 46060201 */  sub.s   $f8, $f0, $f6              
/* 002F0 80897C60 46082282 */  mul.s   $f10, $f4, $f8             
/* 002F4 80897C64 46105481 */  sub.s   $f18, $f10, $f16           
/* 002F8 80897C68 E6120028 */  swc1    $f18, 0x0028($s0)          ## 00000028
.L80897C6C:
/* 002FC 80897C6C 92020160 */  lbu     $v0, 0x0160($s0)           ## 00000160
/* 00300 80897C70 304E0002 */  andi    $t6, $v0, 0x0002           ## $t6 = 00000000
/* 00304 80897C74 51C00022 */  beql    $t6, $zero, .L80897D00     
/* 00308 80897C78 860901B4 */  lh      $t1, 0x01B4($s0)           ## 000001B4
/* 0030C 80897C7C 8607008A */  lh      $a3, 0x008A($s0)           ## 0000008A
/* 00310 80897C80 86040032 */  lh      $a0, 0x0032($s0)           ## 00000032
/* 00314 80897C84 304FFFFC */  andi    $t7, $v0, 0xFFFC           ## $t7 = 00000000
/* 00318 80897C88 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 0031C 80897C8C 00E41823 */  subu    $v1, $a3, $a0              
/* 00320 80897C90 00031C00 */  sll     $v1, $v1, 16               
/* 00324 80897C94 00031C03 */  sra     $v1, $v1, 16               
/* 00328 80897C98 2861C001 */  slti    $at, $v1, 0xC001           
/* 0032C 80897C9C 14200007 */  bne     $at, $zero, .L80897CBC     
/* 00330 80897CA0 A20F0160 */  sb      $t7, 0x0160($s0)           ## 00000160
/* 00334 80897CA4 28614000 */  slti    $at, $v1, 0x4000           
/* 00338 80897CA8 10200004 */  beq     $at, $zero, .L80897CBC     
/* 0033C 80897CAC 34018000 */  ori     $at, $zero, 0x8000         ## $at = 00008000
/* 00340 80897CB0 0081C021 */  addu    $t8, $a0, $at              
/* 00344 80897CB4 A6180032 */  sh      $t8, 0x0032($s0)           ## 00000032
/* 00348 80897CB8 8607008A */  lh      $a3, 0x008A($s0)           ## 0000008A
.L80897CBC:
/* 0034C 80897CBC 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 00350 80897CC0 E7A60010 */  swc1    $f6, 0x0010($sp)           
/* 00354 80897CC4 0C00BDB5 */  jal     func_8002F6D4              
/* 00358 80897CC8 AFA00014 */  sw      $zero, 0x0014($sp)         
/* 0035C 80897CCC 8FB9002C */  lw      $t9, 0x002C($sp)           
/* 00360 80897CD0 2405083E */  addiu   $a1, $zero, 0x083E         ## $a1 = 0000083E
/* 00364 80897CD4 0C00BDF7 */  jal     func_8002F7DC              
/* 00368 80897CD8 8F241C44 */  lw      $a0, 0x1C44($t9)           ## 00001C44
/* 0036C 80897CDC 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 00370 80897CE0 44812000 */  mtc1    $at, $f4                   ## $f4 = 10.00
/* 00374 80897CE4 3C013F00 */  lui     $at, 0x3F00                ## $at = 3F000000
/* 00378 80897CE8 44814000 */  mtc1    $at, $f8                   ## $f8 = 0.50
/* 0037C 80897CEC 24080001 */  addiu   $t0, $zero, 0x0001         ## $t0 = 00000001
/* 00380 80897CF0 A60801B4 */  sh      $t0, 0x01B4($s0)           ## 000001B4
/* 00384 80897CF4 E60401B8 */  swc1    $f4, 0x01B8($s0)           ## 000001B8
/* 00388 80897CF8 E60801B0 */  swc1    $f8, 0x01B0($s0)           ## 000001B0
/* 0038C 80897CFC 860901B4 */  lh      $t1, 0x01B4($s0)           ## 000001B4
.L80897D00:
/* 00390 80897D00 260401B0 */  addiu   $a0, $s0, 0x01B0           ## $a0 = 000001B0
/* 00394 80897D04 3C053F80 */  lui     $a1, 0x3F80                ## $a1 = 3F800000
/* 00398 80897D08 11200019 */  beq     $t1, $zero, .L80897D70     
/* 0039C 80897D0C 3C063D23 */  lui     $a2, 0x3D23                ## $a2 = 3D230000
/* 003A0 80897D10 3C018089 */  lui     $at, %hi(D_80897FE0)       ## $at = 80890000
/* 003A4 80897D14 C4207FE0 */  lwc1    $f0, %lo(D_80897FE0)($at)  
/* 003A8 80897D18 3C013FC0 */  lui     $at, 0x3FC0                ## $at = 3FC00000
/* 003AC 80897D1C 44818000 */  mtc1    $at, $f16                  ## $f16 = 1.50
/* 003B0 80897D20 C60A01B8 */  lwc1    $f10, 0x01B8($s0)          ## 000001B8
/* 003B4 80897D24 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 003B8 80897D28 44812000 */  mtc1    $at, $f4                   ## $f4 = 10.00
/* 003BC 80897D2C 46105481 */  sub.s   $f18, $f10, $f16           
/* 003C0 80897D30 C60A00BC */  lwc1    $f10, 0x00BC($s0)          ## 000000BC
/* 003C4 80897D34 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 003C8 80897D38 E61201B8 */  swc1    $f18, 0x01B8($s0)          ## 000001B8
/* 003CC 80897D3C C60601B8 */  lwc1    $f6, 0x01B8($s0)           ## 000001B8
/* 003D0 80897D40 46043202 */  mul.s   $f8, $f6, $f4              
/* 003D4 80897D44 46085400 */  add.s   $f16, $f10, $f8            
/* 003D8 80897D48 E61000BC */  swc1    $f16, 0x00BC($s0)          ## 000000BC
/* 003DC 80897D4C C61200BC */  lwc1    $f18, 0x00BC($s0)          ## 000000BC
/* 003E0 80897D50 4600903C */  c.lt.s  $f18, $f0                  
/* 003E4 80897D54 00000000 */  nop
/* 003E8 80897D58 45000007 */  bc1f    .L80897D78                 
/* 003EC 80897D5C 00000000 */  nop
/* 003F0 80897D60 0C225F77 */  jal     func_80897DDC              
/* 003F4 80897D64 E60000BC */  swc1    $f0, 0x00BC($s0)           ## 000000BC
/* 003F8 80897D68 10000003 */  beq     $zero, $zero, .L80897D78   
/* 003FC 80897D6C 00000000 */  nop
.L80897D70:
/* 00400 80897D70 0C01DE80 */  jal     Math_ApproxF
              
/* 00404 80897D74 34C6D70A */  ori     $a2, $a2, 0xD70A           ## $a2 = 0000D70A
.L80897D78:
/* 00408 80897D78 3C018089 */  lui     $at, %hi(D_80897FE4)       ## $at = 80890000
/* 0040C 80897D7C C4267FE4 */  lwc1    $f6, %lo(D_80897FE4)($at)  
/* 00410 80897D80 C6020024 */  lwc1    $f2, 0x0024($s0)           ## 00000024
/* 00414 80897D84 240AC000 */  addiu   $t2, $zero, 0xC000         ## $t2 = FFFFC000
/* 00418 80897D88 3C018089 */  lui     $at, %hi(D_80897FE8)       ## $at = 80890000
/* 0041C 80897D8C 4602303C */  c.lt.s  $f6, $f2                   
/* 00420 80897D90 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00424 80897D94 45000003 */  bc1f    .L80897DA4                 
/* 00428 80897D98 00000000 */  nop
/* 0042C 80897D9C 10000008 */  beq     $zero, $zero, .L80897DC0   
/* 00430 80897DA0 A60A0032 */  sh      $t2, 0x0032($s0)           ## 00000032
.L80897DA4:
/* 00434 80897DA4 C4247FE8 */  lwc1    $f4, %lo(D_80897FE8)($at)  
/* 00438 80897DA8 240B4000 */  addiu   $t3, $zero, 0x4000         ## $t3 = 00004000
/* 0043C 80897DAC 4604103C */  c.lt.s  $f2, $f4                   
/* 00440 80897DB0 00000000 */  nop
/* 00444 80897DB4 45000002 */  bc1f    .L80897DC0                 
/* 00448 80897DB8 00000000 */  nop
/* 0044C 80897DBC A60B0032 */  sh      $t3, 0x0032($s0)           ## 00000032
.L80897DC0:
/* 00450 80897DC0 0C00BE0A */  jal     Audio_PlayActorSound2
              
/* 00454 80897DC4 240520B8 */  addiu   $a1, $zero, 0x20B8         ## $a1 = 000020B8
/* 00458 80897DC8 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 0045C 80897DCC 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00460 80897DD0 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 00464 80897DD4 03E00008 */  jr      $ra                        
/* 00468 80897DD8 00000000 */  nop


