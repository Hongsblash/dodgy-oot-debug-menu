glabel func_80B86AC8
/* 001A8 80B86AC8 27BDFFD8 */  addiu   $sp, $sp, 0xFFD8           ## $sp = FFFFFFD8
/* 001AC 80B86ACC AFBF0024 */  sw      $ra, 0x0024($sp)           
/* 001B0 80B86AD0 AFB00020 */  sw      $s0, 0x0020($sp)           
/* 001B4 80B86AD4 00808025 */  or      $s0, $a0, $zero            ## $s0 = 00000000
/* 001B8 80B86AD8 0C00B638 */  jal     Actor_MoveForward
              
/* 001BC 80B86ADC AFA5002C */  sw      $a1, 0x002C($sp)           
/* 001C0 80B86AE0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 001C4 80B86AE4 0C00BD04 */  jal     func_8002F410              
/* 001C8 80B86AE8 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 001CC 80B86AEC 10400005 */  beq     $v0, $zero, .L80B86B04     
/* 001D0 80B86AF0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 001D4 80B86AF4 0C00B55C */  jal     Actor_Kill
              
/* 001D8 80B86AF8 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 001DC 80B86AFC 1000002E */  beq     $zero, $zero, .L80B86BB8   
/* 001E0 80B86B00 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80B86B04:
/* 001E4 80B86B04 3C014248 */  lui     $at, 0x4248                ## $at = 42480000
/* 001E8 80B86B08 44812000 */  mtc1    $at, $f4                   ## $f4 = 50.00
/* 001EC 80B86B0C 8FA5002C */  lw      $a1, 0x002C($sp)           
/* 001F0 80B86B10 24060029 */  addiu   $a2, $zero, 0x0029         ## $a2 = 00000029
/* 001F4 80B86B14 3C0741F0 */  lui     $a3, 0x41F0                ## $a3 = 41F00000
/* 001F8 80B86B18 0C00BD0D */  jal     func_8002F434              
/* 001FC 80B86B1C E7A40010 */  swc1    $f4, 0x0010($sp)           
/* 00200 80B86B20 3C014120 */  lui     $at, 0x4120                ## $at = 41200000
/* 00204 80B86B24 44810000 */  mtc1    $at, $f0                   ## $f0 = 10.00
/* 00208 80B86B28 44803000 */  mtc1    $zero, $f6                 ## $f6 = 0.00
/* 0020C 80B86B2C 240E0005 */  addiu   $t6, $zero, 0x0005         ## $t6 = 00000005
/* 00210 80B86B30 44060000 */  mfc1    $a2, $f0                   
/* 00214 80B86B34 44070000 */  mfc1    $a3, $f0                   
/* 00218 80B86B38 AFAE0014 */  sw      $t6, 0x0014($sp)           
/* 0021C 80B86B3C 8FA4002C */  lw      $a0, 0x002C($sp)           
/* 00220 80B86B40 02002825 */  or      $a1, $s0, $zero            ## $a1 = 00000000
/* 00224 80B86B44 0C00B92D */  jal     func_8002E4B4              
/* 00228 80B86B48 E7A60010 */  swc1    $f6, 0x0010($sp)           
/* 0022C 80B86B4C 960F0088 */  lhu     $t7, 0x0088($s0)           ## 00000088
/* 00230 80B86B50 31F80001 */  andi    $t8, $t7, 0x0001           ## $t8 = 00000000
/* 00234 80B86B54 53000018 */  beql    $t8, $zero, .L80B86BB8     
/* 00238 80B86B58 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 0023C 80B86B5C 8619019A */  lh      $t9, 0x019A($s0)           ## 0000019A
/* 00240 80B86B60 2728FFFF */  addiu   $t0, $t9, 0xFFFF           ## $t0 = FFFFFFFF
/* 00244 80B86B64 A608019A */  sh      $t0, 0x019A($s0)           ## 0000019A
/* 00248 80B86B68 8602019A */  lh      $v0, 0x019A($s0)           ## 0000019A
/* 0024C 80B86B6C 2841003C */  slti    $at, $v0, 0x003C           
/* 00250 80B86B70 1020000C */  beq     $at, $zero, .L80B86BA4     
/* 00254 80B86B74 30490001 */  andi    $t1, $v0, 0x0001           ## $t1 = 00000000
/* 00258 80B86B78 51200007 */  beql    $t1, $zero, .L80B86B98     
/* 0025C 80B86B7C 860C019C */  lh      $t4, 0x019C($s0)           ## 0000019C
/* 00260 80B86B80 860A019C */  lh      $t2, 0x019C($s0)           ## 0000019C
/* 00264 80B86B84 8602019A */  lh      $v0, 0x019A($s0)           ## 0000019A
/* 00268 80B86B88 354B0002 */  ori     $t3, $t2, 0x0002           ## $t3 = 00000002
/* 0026C 80B86B8C 10000005 */  beq     $zero, $zero, .L80B86BA4   
/* 00270 80B86B90 A60B019C */  sh      $t3, 0x019C($s0)           ## 0000019C
/* 00274 80B86B94 860C019C */  lh      $t4, 0x019C($s0)           ## 0000019C
.L80B86B98:
/* 00278 80B86B98 8602019A */  lh      $v0, 0x019A($s0)           ## 0000019A
/* 0027C 80B86B9C 318DFFFD */  andi    $t5, $t4, 0xFFFD           ## $t5 = 00000000
/* 00280 80B86BA0 A60D019C */  sh      $t5, 0x019C($s0)           ## 0000019C
.L80B86BA4:
/* 00284 80B86BA4 54400004 */  bnel    $v0, $zero, .L80B86BB8     
/* 00288 80B86BA8 8FBF0024 */  lw      $ra, 0x0024($sp)           
/* 0028C 80B86BAC 0C00B55C */  jal     Actor_Kill
              
/* 00290 80B86BB0 02002025 */  or      $a0, $s0, $zero            ## $a0 = 00000000
/* 00294 80B86BB4 8FBF0024 */  lw      $ra, 0x0024($sp)           
.L80B86BB8:
/* 00298 80B86BB8 8FB00020 */  lw      $s0, 0x0020($sp)           
/* 0029C 80B86BBC 27BD0028 */  addiu   $sp, $sp, 0x0028           ## $sp = 00000000
/* 002A0 80B86BC0 03E00008 */  jr      $ra                        
/* 002A4 80B86BC4 00000000 */  nop


