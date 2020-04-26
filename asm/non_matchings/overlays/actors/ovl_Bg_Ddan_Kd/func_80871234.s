.rdata
glabel D_80871920
    .asciz "dam    %d\n"
    .balign 4

.text
glabel func_80871234
/* 00144 80871234 27BDFFD0 */  addiu   $sp, $sp, 0xFFD0           ## $sp = FFFFFFD0
/* 00148 80871238 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 0014C 8087123C 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 00150 80871240 AFA50034 */  sw      $a1, 0x0034($sp)           
/* 00154 80871244 00A02025 */  or      $a0, $a1, $zero            ## $a0 = 00000000
/* 00158 80871248 AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 0015C 8087124C 26050178 */  addiu   $a1, $s0, 0x0178           ## $a1 = 00000178
/* 00160 80871250 0C00CD90 */  jal     func_80033640              
/* 00164 80871254 AFA50028 */  sw      $a1, 0x0028($sp)           
/* 00168 80871258 10400009 */  beq     $v0, $zero, .L80871280     
/* 0016C 8087125C 00401825 */  or      $v1, $v0, $zero            ## $v1 = 00000000
/* 00170 80871260 3C048087 */  lui     $a0, %hi(D_80871920)       ## $a0 = 80870000
/* 00174 80871264 920500B0 */  lbu     $a1, 0x00B0($s0)           ## 000000B0
/* 00178 80871268 AFA2002C */  sw      $v0, 0x002C($sp)           
/* 0017C 8087126C 0C00084C */  jal     osSyncPrintf
              
/* 00180 80871270 24841920 */  addiu   $a0, $a0, %lo(D_80871920)  ## $a0 = 80871920
/* 00184 80871274 8FA3002C */  lw      $v1, 0x002C($sp)           
/* 00188 80871278 240E0002 */  addiu   $t6, $zero, 0x0002         ## $t6 = 00000002
/* 0018C 8087127C A46E001C */  sh      $t6, 0x001C($v1)           ## 0000001C
.L80871280:
/* 00190 80871280 5060001C */  beql    $v1, $zero, .L808712F4     
/* 00194 80871284 86020168 */  lh      $v0, 0x0168($s0)           ## 00000168
/* 00198 80871288 8E020164 */  lw      $v0, 0x0164($s0)           ## 00000164
/* 0019C 8087128C 50400019 */  beql    $v0, $zero, .L808712F4     
/* 001A0 80871290 86020168 */  lh      $v0, 0x0168($s0)           ## 00000168
/* 001A4 80871294 10620016 */  beq     $v1, $v0, .L808712F0       
/* 001A8 80871298 2604016C */  addiu   $a0, $s0, 0x016C           ## $a0 = 0000016C
/* 001AC 8087129C 24650024 */  addiu   $a1, $v1, 0x0024           ## $a1 = 00000024
/* 001B0 808712A0 0C01E00A */  jal     Math_Vec3f_DistXZ
              
/* 001B4 808712A4 AFA3002C */  sw      $v1, 0x002C($sp)           
/* 001B8 808712A8 3C0142A0 */  lui     $at, 0x42A0                ## $at = 42A00000
/* 001BC 808712AC 44812000 */  mtc1    $at, $f4                   ## $f4 = 80.00
/* 001C0 808712B0 8FA3002C */  lw      $v1, 0x002C($sp)           
/* 001C4 808712B4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 001C8 808712B8 4600203C */  c.lt.s  $f4, $f0                   
/* 001CC 808712BC 3C058087 */  lui     $a1, %hi(func_80871364)    ## $a1 = 80870000
/* 001D0 808712C0 4502000C */  bc1fl   .L808712F4                 
/* 001D4 808712C4 86020168 */  lh      $v0, 0x0168($s0)           ## 00000168
/* 001D8 808712C8 0C21C43C */  jal     BgDdanKd_SetupAction              
/* 001DC 808712CC 24A51364 */  addiu   $a1, $a1, %lo(func_80871364) ## $a1 = 80871364
/* 001E0 808712D0 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 001E4 808712D4 24050BEA */  addiu   $a1, $zero, 0x0BEA         ## $a1 = 00000BEA
/* 001E8 808712D8 240603E7 */  addiu   $a2, $zero, 0x03E7         ## $a2 = 000003E7
/* 001EC 808712DC 02003825 */  or      $a3, $s0, $zero            ## $a3 = 00000000
/* 001F0 808712E0 0C02003E */  jal     func_800800F8              
/* 001F4 808712E4 AFA00010 */  sw      $zero, 0x0010($sp)         
/* 001F8 808712E8 1000001A */  beq     $zero, $zero, .L80871354   
/* 001FC 808712EC 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L808712F0:
/* 00200 808712F0 86020168 */  lh      $v0, 0x0168($s0)           ## 00000168
.L808712F4:
/* 00204 808712F4 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00208 808712F8 10400003 */  beq     $v0, $zero, .L80871308     
/* 0020C 808712FC 244FFFFF */  addiu   $t7, $v0, 0xFFFF           ## $t7 = FFFFFFFF
/* 00210 80871300 1000000B */  beq     $zero, $zero, .L80871330   
/* 00214 80871304 A60F0168 */  sh      $t7, 0x0168($s0)           ## 00000168
.L80871308:
/* 00218 80871308 10600009 */  beq     $v1, $zero, .L80871330     
/* 0021C 8087130C AE030164 */  sw      $v1, 0x0164($s0)           ## 00000164
/* 00220 80871310 2418000D */  addiu   $t8, $zero, 0x000D         ## $t8 = 0000000D
/* 00224 80871314 A6180168 */  sh      $t8, 0x0168($s0)           ## 00000168
/* 00228 80871318 8C680024 */  lw      $t0, 0x0024($v1)           ## 00000024
/* 0022C 8087131C AE08016C */  sw      $t0, 0x016C($s0)           ## 0000016C
/* 00230 80871320 8C790028 */  lw      $t9, 0x0028($v1)           ## 00000028
/* 00234 80871324 AE190170 */  sw      $t9, 0x0170($s0)           ## 00000170
/* 00238 80871328 8C68002C */  lw      $t0, 0x002C($v1)           ## 0000002C
/* 0023C 8087132C AE080174 */  sw      $t0, 0x0174($s0)           ## 00000174
.L80871330:
/* 00240 80871330 0C0189B7 */  jal     Collider_CylinderUpdate
              
/* 00244 80871334 8FA50028 */  lw      $a1, 0x0028($sp)           
/* 00248 80871338 8FA40034 */  lw      $a0, 0x0034($sp)           
/* 0024C 8087133C 3C010001 */  lui     $at, 0x0001                ## $at = 00010000
/* 00250 80871340 34211E60 */  ori     $at, $at, 0x1E60           ## $at = 00011E60
/* 00254 80871344 8FA60028 */  lw      $a2, 0x0028($sp)           
/* 00258 80871348 0C01767D */  jal     CollisionCheck_SetAC
              ## CollisionCheck_setAC
/* 0025C 8087134C 00812821 */  addu    $a1, $a0, $at              
/* 00260 80871350 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80871354:
/* 00264 80871354 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 00268 80871358 27BD0030 */  addiu   $sp, $sp, 0x0030           ## $sp = 00000000
/* 0026C 8087135C 03E00008 */  jr      $ra                        
/* 00270 80871360 00000000 */  nop
