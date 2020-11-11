clear all
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook:
%    C:\Users\laura.mihai\Desktop\SENSECO-Training\15ms\5fL_15ms.xlsx
%    Worksheet: 1801064U1
%
% To extend the code for use with different selected data or a different
% spreadsheet, generate a function instead of a script.

% Auto-generated by MATLAB on 2019/10/27 12:16:37

%% write wavelength values
Wavelength=[623.84;623.97;624.09;624.21;624.34;624.46;624.59;624.71;624.84;624.96;625.08;625.21;625.33;625.46;625.58;625.71;625.83;625.95;626.08;626.2;626.33;626.45;626.57;626.7;626.82;626.95;627.07;627.19;627.32;627.44;627.57;627.69;627.81;627.94;628.06;628.19;628.31;628.43;628.56;628.68;628.8;628.93;629.05;629.17;629.3;629.42;629.55;629.67;629.79;629.92;630.04;630.16;630.29;630.41;630.53;630.66;630.78;630.9;631.03;631.15;631.27;631.4;631.52;631.64;631.77;631.89;632.01;632.14;632.26;632.38;632.51;632.63;632.75;632.88;633;633.12;633.25;633.37;633.49;633.62;633.74;633.86;633.98;634.11;634.23;634.35;634.48;634.6;634.72;634.85;634.97;635.09;635.21;635.34;635.46;635.58;635.71;635.83;635.95;636.07;636.2;636.32;636.44;636.56;636.69;636.81;636.93;637.05;637.18;637.3;637.42;637.54;637.67;637.79;637.91;638.03;638.16;638.28;638.4;638.52;638.65;638.77;638.89;639.01;639.14;639.26;639.38;639.5;639.62;639.75;639.87;639.99;640.11;640.24;640.36;640.48;640.6;640.72;640.85;640.97;641.09;641.21;641.33;641.46;641.58;641.7;641.82;641.94;642.07;642.19;642.31;642.43;642.55;642.67;642.8;642.92;643.04;643.16;643.28;643.41;643.53;643.65;643.77;643.89;644.01;644.14;644.26;644.38;644.5;644.62;644.74;644.86;644.99;645.11;645.23;645.35;645.47;645.59;645.71;645.84;645.96;646.08;646.2;646.32;646.44;646.56;646.68;646.81;646.93;647.05;647.17;647.29;647.41;647.53;647.65;647.77;647.9;648.02;648.14;648.26;648.38;648.5;648.62;648.74;648.86;648.98;649.11;649.23;649.35;649.47;649.59;649.71;649.83;649.95;650.07;650.19;650.31;650.43;650.55;650.67;650.8;650.92;651.04;651.16;651.28;651.4;651.52;651.64;651.76;651.88;652;652.12;652.24;652.36;652.48;652.6;652.72;652.84;652.96;653.08;653.2;653.32;653.44;653.57;653.69;653.81;653.93;654.05;654.17;654.29;654.41;654.53;654.65;654.77;654.89;655.01;655.13;655.25;655.37;655.49;655.61;655.73;655.85;655.97;656.09;656.21;656.33;656.45;656.56;656.68;656.8;656.92;657.04;657.16;657.28;657.4;657.52;657.64;657.76;657.88;658;658.12;658.24;658.36;658.48;658.6;658.72;658.84;658.96;659.08;659.2;659.31;659.43;659.55;659.67;659.79;659.91;660.03;660.15;660.27;660.39;660.51;660.63;660.75;660.87;660.98;661.1;661.22;661.34;661.46;661.58;661.7;661.82;661.94;662.06;662.17;662.29;662.41;662.53;662.65;662.77;662.89;663.01;663.13;663.24;663.36;663.48;663.6;663.72;663.84;663.96;664.08;664.19;664.31;664.43;664.55;664.67;664.79;664.91;665.02;665.14;665.26;665.38;665.5;665.62;665.73;665.85;665.97;666.09;666.21;666.33;666.44;666.56;666.68;666.8;666.92;667.04;667.15;667.27;667.39;667.51;667.63;667.74;667.86;667.98;668.1;668.22;668.34;668.45;668.57;668.69;668.81;668.92;669.04;669.16;669.28;669.4;669.51;669.63;669.75;669.87;669.99;670.1;670.22;670.34;670.46;670.57;670.69;670.81;670.93;671.04;671.16;671.28;671.4;671.51;671.63;671.75;671.87;671.98;672.1;672.22;672.34;672.45;672.57;672.69;672.81;672.92;673.04;673.16;673.28;673.39;673.51;673.63;673.74;673.86;673.98;674.1;674.21;674.33;674.45;674.56;674.68;674.8;674.92;675.03;675.15;675.27;675.38;675.5;675.62;675.73;675.85;675.97;676.08;676.2;676.32;676.43;676.55;676.67;676.78;676.9;677.02;677.13;677.25;677.37;677.48;677.6;677.72;677.83;677.95;678.07;678.18;678.3;678.42;678.53;678.65;678.77;678.88;679;679.11;679.23;679.35;679.46;679.58;679.7;679.81;679.93;680.04;680.16;680.28;680.39;680.51;680.63;680.74;680.86;680.97;681.09;681.21;681.32;681.44;681.55;681.67;681.79;681.9;682.02;682.13;682.25;682.36;682.48;682.6;682.71;682.83;682.94;683.06;683.18;683.29;683.41;683.52;683.64;683.75;683.87;683.98;684.1;684.22;684.33;684.45;684.56;684.68;684.79;684.91;685.02;685.14;685.25;685.37;685.48;685.6;685.72;685.83;685.95;686.06;686.18;686.29;686.41;686.52;686.64;686.75;686.87;686.98;687.1;687.21;687.33;687.44;687.56;687.67;687.79;687.9;688.02;688.13;688.25;688.36;688.48;688.59;688.71;688.82;688.94;689.05;689.17;689.28;689.39;689.51;689.62;689.74;689.85;689.97;690.08;690.2;690.31;690.43;690.54;690.66;690.77;690.88;691;691.11;691.23;691.34;691.46;691.57;691.68;691.8;691.91;692.03;692.14;692.26;692.37;692.48;692.6;692.71;692.83;692.94;693.06;693.17;693.28;693.4;693.51;693.63;693.74;693.85;693.97;694.08;694.19;694.31;694.42;694.54;694.65;694.76;694.88;694.99;695.11;695.22;695.33;695.45;695.56;695.67;695.79;695.9;696.01;696.13;696.24;696.36;696.47;696.58;696.7;696.81;696.92;697.04;697.15;697.26;697.38;697.49;697.6;697.72;697.83;697.94;698.06;698.17;698.28;698.4;698.51;698.62;698.74;698.85;698.96;699.07;699.19;699.3;699.41;699.53;699.64;699.75;699.87;699.98;700.09;700.2;700.32;700.43;700.54;700.66;700.77;700.88;700.99;701.11;701.22;701.33;701.45;701.56;701.67;701.78;701.9;702.01;702.12;702.23;702.35;702.46;702.57;702.68;702.8;702.91;703.02;703.13;703.25;703.36;703.47;703.58;703.69;703.81;703.92;704.03;704.14;704.26;704.37;704.48;704.59;704.7;704.82;704.93;705.04;705.15;705.27;705.38;705.49;705.6;705.71;705.82;705.94;706.05;706.16;706.27;706.38;706.5;706.61;706.72;706.83;706.94;707.06;707.17;707.28;707.39;707.5;707.61;707.72;707.84;707.95;708.06;708.17;708.28;708.39;708.51;708.62;708.73;708.84;708.95;709.06;709.17;709.29;709.4;709.51;709.62;709.73;709.84;709.95;710.06;710.18;710.29;710.4;710.51;710.62;710.73;710.84;710.95;711.06;711.18;711.29;711.4;711.51;711.62;711.73;711.84;711.95;712.06;712.17;712.28;712.4;712.51;712.62;712.73;712.84;712.95;713.06;713.17;713.28;713.39;713.5;713.61;713.72;713.83;713.94;714.05;714.17;714.28;714.39;714.5;714.61;714.72;714.83;714.94;715.05;715.16;715.27;715.38;715.49;715.6;715.71;715.82;715.93;716.04;716.15;716.26;716.37;716.48;716.59;716.7;716.81;716.92;717.03;717.14;717.25;717.36;717.47;717.58;717.69;717.8;717.91;718.02;718.13;718.24;718.35;718.46;718.57;718.68;718.79;718.9;719.01;719.12;719.23;719.34;719.45;719.56;719.66;719.77;719.88;719.99;720.1;720.21;720.32;720.43;720.54;720.65;720.76;720.87;720.98;721.09;721.2;721.31;721.41;721.52;721.63;721.74;721.85;721.96;722.07;722.18;722.29;722.4;722.51;722.61;722.72;722.83;722.94;723.05;723.16;723.27;723.38;723.49;723.59;723.7;723.81;723.92;724.03;724.14;724.25;724.36;724.46;724.57;724.68;724.79;724.9;725.01;725.12;725.22;725.33;725.44;725.55;725.66;725.77;725.87;725.98;726.09;726.2;726.31;726.42;726.52;726.63;726.74;726.85;726.96;727.07;727.17;727.28;727.39;727.5;727.61;727.71;727.82;727.93;728.04;728.15;728.25;728.36;728.47;728.58;728.69;728.79;728.9;729.01;729.12;729.22;729.33;729.44;729.55;729.66;729.76;729.87;729.98;730.09;730.19;730.3;730.41;730.52;730.62;730.73;730.84;730.95;731.05;731.16;731.27;731.38;731.48;731.59;731.7;731.8;731.91;732.02;732.13;732.23;732.34;732.45;732.56;732.66;732.77;732.88;732.98;733.09;733.2;733.3;733.41;733.52;733.63;733.73;733.84;733.95;734.05;734.16;734.27;734.37;734.48;734.59;734.69;734.8;734.91;735.01;735.12;735.23;735.33;735.44;735.55;735.65;735.76;735.87;735.97;736.08;736.19;736.29;736.4;736.51;736.61;736.72;736.82;736.93;737.04;737.14;737.25;737.36;737.46;737.57;737.67;737.78;737.89;737.99;738.1;738.2;738.31;738.42;738.52;738.63;738.73;738.84;738.95;739.05;739.16;739.26;739.37;739.48;739.58;739.69;739.79;739.9;740;740.11;740.22;740.32;740.43;740.53;740.64;740.74;740.85;740.96;741.06;741.17;741.27;741.38;741.48;741.59;741.69;741.8;741.9;742.01;742.12;742.22;742.33;742.43;742.54;742.64;742.75;742.85;742.96;743.06;743.17;743.27;743.38;743.48;743.59;743.69;743.8;743.9;744.01;744.11;744.22;744.32;744.43;744.53;744.64;744.74;744.85;744.95;745.06;745.16;745.27;745.37;745.47;745.58;745.68;745.79;745.89;746;746.1;746.21;746.31;746.42;746.52;746.62;746.73;746.83;746.94;747.04;747.15;747.25;747.36;747.46;747.56;747.67;747.77;747.88;747.98;748.08;748.19;748.29;748.4;748.5;748.6;748.71;748.81;748.92;749.02;749.12;749.23;749.33;749.44;749.54;749.64;749.75;749.85;749.96;750.06;750.16;750.27;750.37;750.47;750.58;750.68;750.79;750.89;750.99;751.1;751.2;751.3;751.41;751.51;751.61;751.72;751.82;751.92;752.03;752.13;752.23;752.34;752.44;752.54;752.65;752.75;752.85;752.96;753.06;753.16;753.27;753.37;753.47;753.57;753.68;753.78;753.88;753.99;754.09;754.19;754.3;754.4;754.5;754.6;754.71;754.81;754.91;755.02;755.12;755.22;755.32;755.43;755.53;755.63;755.73;755.84;755.94;756.04;756.14;756.25;756.35;756.45;756.55;756.66;756.76;756.86;756.96;757.07;757.17;757.27;757.37;757.48;757.58;757.68;757.78;757.88;757.99;758.09;758.19;758.29;758.39;758.5;758.6;758.7;758.8;758.9;759.01;759.11;759.21;759.31;759.41;759.52;759.62;759.72;759.82;759.92;760.02;760.13;760.23;760.33;760.43;760.53;760.63;760.74;760.84;760.94;761.04;761.14;761.24;761.35;761.45;761.55;761.65;761.75;761.85;761.95;762.06;762.16;762.26;762.36;762.46;762.56;762.66;762.76;762.86;762.97;763.07;763.17;763.27;763.37;763.47;763.57;763.67;763.77;763.88;763.98;764.08;764.18;764.28;764.38;764.48;764.58;764.68;764.78;764.88;764.98;765.08;765.19;765.29;765.39;765.49;765.59;765.69;765.79;765.89;765.99;766.09;766.19;766.29;766.39;766.49;766.59;766.69;766.79;766.89;766.99;767.09;767.19;767.29;767.39;767.49;767.59;767.69;767.79;767.89;767.99;768.09;768.2;768.3;768.39;768.49;768.59;768.69;768.79;768.89;768.99;769.09;769.19;769.29;769.39;769.49;769.59;769.69;769.79;769.89;769.99;770.09;770.19;770.29;770.39;770.49;770.59;770.69;770.79;770.89;770.99;771.09;771.19;771.28;771.38;771.48;771.58;771.68;771.78;771.88;771.98;772.08;772.18;772.28;772.38;772.48;772.57;772.67;772.77;772.87;772.97;773.07;773.17;773.27;773.37;773.47;773.56;773.66;773.76;773.86;773.96;774.06;774.16;774.26;774.35;774.45;774.55;774.65;774.75;774.85;774.95;775.04;775.14;775.24;775.34;775.44;775.54;775.64;775.73;775.83;775.93;776.03;776.13;776.23;776.32;776.42;776.52;776.62;776.72;776.81;776.91;777.01;777.11;777.21;777.3;777.4;777.5;777.6;777.7;777.79;777.89;777.99;778.09;778.19;778.28;778.38;778.48;778.58;778.67;778.77;778.87;778.97;779.07;779.16;779.26;779.36;779.46;779.55;779.65;779.75;779.85;779.94;780.04;780.14;780.24;780.33;780.43;780.53;780.63;780.72;780.82;780.92;781.01;781.11;781.21;781.31;781.4;781.5;781.6;781.69;781.79;781.89;781.98;782.08;782.18;782.28;782.37;782.47;782.57;782.66;782.76;782.86;782.95;783.05;783.15;783.24;783.34;783.44;783.53;783.63;783.73;783.82;783.92;784.02;784.11;784.21;784.31;784.4;784.5;784.6;784.69;784.79;784.88;784.98;785.08;785.17;785.27;785.37;785.46;785.56;785.65;785.75;785.85;785.94;786.04;786.13;786.23;786.33;786.42;786.52;786.61;786.71;786.81;786.9;787;787.09;787.19;787.29;787.38;787.48;787.57;787.67;787.76;787.86;787.96;788.05;788.15;788.24;788.34;788.43;788.53;788.62;788.72;788.82;788.91;789.01;789.1;789.2;789.29;789.39;789.48;789.58;789.67;789.77;789.86;789.96;790.05;790.15;790.24;790.34;790.43;790.53;790.62;790.72;790.81;790.91;791;791.1;791.19;791.29;791.38;791.48;791.57;791.67;791.76;791.86;791.95;792.05;792.14;792.24;792.33;792.43;792.52;792.61;792.71;792.8;792.9;792.99;793.09;793.18;793.28;793.37;793.46;793.56;793.65;793.75;793.84;793.94;794.03;794.12;794.22;794.31;794.41;794.5;794.59;794.69;794.78;794.88;794.97;795.06;795.16;795.25;795.35;795.44;795.53;795.63;795.72;795.82;795.91;796;796.1;796.19;796.28;796.38;796.47;796.57;796.66;796.75;796.85;796.94;797.03;797.13;797.22;797.31;797.41;797.5;797.59;797.69;797.78;797.87;797.97;798.06;798.15;798.25;798.34;798.43;798.53;798.62;798.71;798.81;798.9;798.99;799.08;799.18;799.27;799.36;799.46;799.55;799.64;799.73;799.83;799.92;800.01;800.11;800.2;800.29;800.38;800.48;800.57;800.66;800.75;800.85;800.94;801.03;801.12;801.22;801.31;801.4;801.49;801.59;801.68;801.77;801.86;801.96;802.05;802.14;802.23;802.33;802.42;802.51;802.6;802.69;802.79;802.88;802.97;803.06;803.15;803.25;803.34;803.43;803.52;803.61;803.71;803.8;803.89;803.98;804.07;804.16;804.26;804.35;804.44;804.53;804.62;804.71;804.81;804.9;804.99;805.08;805.17;805.26;805.36;805.45;805.54;805.63;805.72;805.81;805.9;806;806.09;806.18;806.27;806.36;806.45;806.54;806.63;806.72;806.82;806.91;807;807.09;807.18;807.27;807.36;807.45;807.54;807.63;807.73;807.82;807.91;808;808.09;808.18;808.27;808.36;808.45;808.54;808.63;808.72;808.81;808.9;809;809.09;809.18;809.27;809.36;809.45;809.54;809.63;809.72;809.81;809.9;809.99;810.08;810.17;810.26;810.35;810.44;810.53;810.62;810.71;810.8;810.89;810.98;811.07;811.16;811.25;811.34;811.43;811.52;811.61;811.7;811.79;811.88;811.97;812.06;812.15;812.24;812.33;812.42;812.51;812.6;812.69;812.78;812.87;812.96;813.05;813.14;813.23;813.32;813.4;813.49;813.58;813.67;813.76;813.85;813.94;814.03;814.12;814.21;814.3;814.39;814.48;814.56;814.65;814.74;814.83;814.92;815.01;815.1;815.19;815.28;815.37;815.45;815.54;815.63;815.72;815.81;815.9;815.99;816.08;816.17;816.25;816.34;816.43;816.52;816.61;816.7;816.79;816.87;816.96;817.05;817.14;817.23;817.32;817.4;817.49;817.58;817.67;817.76;817.85;817.93;818.02;818.11;818.2;818.29;818.38;818.46;818.55;818.64;818.73;818.82;818.9;818.99;819.08;819.17;819.26;819.34;819.43;819.52;819.61;819.69;819.78;819.87;819.96;820.05;820.13;820.22;820.31;820.4;820.48;820.57;820.66;820.75;820.83;820.92;821.01;821.1;821.18;821.27;821.36;821.45;821.53;821.62;821.71;821.79;821.88;821.97;822.06;822.14;822.23;822.32;822.4;822.49;822.58;822.67;822.75;822.84;822.93;823.01;823.1;823.19;823.27;823.36;823.45;823.53;823.62;823.71;823.79;823.88;823.97;824.05;824.14;824.23;824.31;824.4;824.49;824.57;824.66;824.75;824.83;824.92;825.01;825.09;825.18;825.26;825.35;825.44;825.52;825.61;825.7;825.78;825.87;825.95;826.04;826.13;826.21;826.3;826.38;826.47;826.56;826.64;826.73;826.81;826.9;826.98;827.07;827.16;827.24;827.33;827.41;827.5;827.58;827.67;827.76;827.84;827.93;828.01;828.1;828.18;828.27;828.35;828.44;828.53;828.61;828.7;828.78;828.87;828.95;829.04;829.12;829.21;829.29;829.38;829.46;829.55;829.63;829.72;829.8;829.89;829.97;830.06;830.14;830.23;830.31;830.4;830.48;830.57;830.65;830.74;830.82;830.91;830.99;831.08;831.16;831.24;831.33;831.41;831.5;831.58;831.67;831.75;831.84;831.92;832.01;832.09;832.17;832.26;832.34;832.43;832.51;832.6;832.68;832.76;832.85;832.93;833.02;833.1;833.18;833.27;833.35;833.44;833.52;833.6;833.69;833.77;833.86;833.94;834.02;834.11;834.19;834.28;834.36;834.44;834.53;834.61;834.69;834.78;834.86;834.94;835.03;835.11;835.2;835.28;835.36;835.45;835.53;835.61;835.7;835.78;835.86;835.95;836.03;836.11;836.2;836.28;836.36;836.45;836.53;836.61;836.69;836.78;836.86;836.94;837.03;837.11;837.19;837.28;837.36;837.44;837.52;837.61;837.69;837.77;837.85;837.94;838.02;838.1;838.19
];

%% Import the data, extracting spreadsheet dates in MATLAB serial date number format (datenum)
%% read light data level 5fL
[~, ~, raw, dateNums] = xlsread('C:\Users\laura.mihai\Desktop\SENSECO-Training\15ms\5fL_15ms.xlsx','1801064U1','B6:DQ2052','',@convertSpreadsheetDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
%% Replace date strings by MATLAB serial date numbers (datenum)
R = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); % Find spreadsheet dates
raw(R) = dateNums(R);
%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
fL_5 = reshape([raw{:}],size(raw));
Lg5=fL_5(1:2047,1:100);
MeanLg5=mean(Lg5,2);
StdevLg5=std(Lg5,0,2);


%% read light data level 100fL
[~, ~, raw, dateNums] = xlsread('C:\Users\laura.mihai\Desktop\SENSECO-Training\15ms\100fL_15ms.xlsx','1801064U1','B6:DQ2052','',@convertSpreadsheetDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
%% Replace date strings by MATLAB serial date numbers (datenum)
R = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); % Find spreadsheet dates
raw(R) = dateNums(R);
%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
fL_100 = reshape([raw{:}],size(raw));
Lg100=fL_100(1:2047,1:100);
MeanLg100=mean(Lg100,2);
StdevLg100=std(Lg100,0,2);

%% read light data level 1000fL
[~, ~, raw, dateNums] = xlsread('C:\Users\laura.mihai\Desktop\SENSECO-Training\15ms\1000fL_15ms.xlsx','1801064U1','B6:DQ2052','',@convertSpreadsheetDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
%% Replace date strings by MATLAB serial date numbers (datenum)
R = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); % Find spreadsheet dates
raw(R) = dateNums(R);
%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
fL_1000 = reshape([raw{:}],size(raw));
Lg1000=fL_1000(1:2047,1:100);
MeanLg1000=mean(Lg1000,2);
StdevLg1000=std(Lg1000,0,2);

%% read light data level 10000fL
[~, ~, raw, dateNums] = xlsread('C:\Users\laura.mihai\Desktop\SENSECO-Training\15ms\10000fL_15ms.xlsx','1801064U1','B6:DQ2052','',@convertSpreadsheetDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
%% Replace date strings by MATLAB serial date numbers (datenum)
R = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); % Find spreadsheet dates
raw(R) = dateNums(R);
%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
fL_10000 = reshape([raw{:}],size(raw));
Lg10000=fL_10000(1:2047,1:100);
MeanLg10000=mean(Lg10000,2);
StdevLg10000=std(Lg10000,0,2);

%% Read dark signal from the begining 
[~, ~, raw, dateNums] = xlsread('C:\Users\laura.mihai\Desktop\SENSECO-Training\15ms\Dark15ms.xlsx','1801064U1','B6:HF2053','',@convertSpreadsheetDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

%% Replace date strings by MATLAB serial date numbers (datenum)
R = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); % Find spreadsheet dates
raw(R) = dateNums(R);

%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
Dark15ms = reshape([raw{:}],size(raw));
Dk=Dark15ms(1:2047,1:100);
MeanDk=mean(Dk,2);
StdevDk=std(Dk,0,2);

DC5=MeanLg5-MeanDk;
DC100=MeanLg100-MeanDk;
DC1000=MeanLg1000-MeanDk;
DC10000=MeanLg10000-MeanDk;

%% plots
figure
subplot(2,2,1);
plot(Wavelength,MeanLg5,'b',Wavelength,MeanLg100,'y',Wavelength,MeanLg1000,'r',Wavelength,MeanLg10000,'c')
title('Light Signal Mean')
xlabel('Wavelength [nm]')
ylabel('DN light')

subplot(2,2,3);
plot(Wavelength,StdevLg5,'b--',Wavelength,StdevLg100,'y--',Wavelength,StdevLg1000,'r--',Wavelength,StdevLg10000,'c--')
title('Light Signal Standard Deviation')
xlabel('Wavelength [nm]')
ylabel('Stdev light')

subplot(2,2,2);
plot(Wavelength,MeanDk,'b')
title('Dark Signal Mean')
xlabel('Wavelength [nm]')
ylabel('DN Dark')

subplot(2,2,4);
plot(Wavelength,StdevDk,'b--')
title('Dark Signal Standard Deviation')
xlabel('Wavelength [nm]')
ylabel('Stdev Dark')


figure
plot(Wavelength,DC5,'b--',Wavelength,DC100,'g--', Wavelength,DC1000,'y--',Wavelength,DC10000,'c--')
title('DC correction')
xlabel('Wavelength [nm]')
ylabel('DC[DN]')
%% Clear temporary variables
clearvars raw dateNums R;