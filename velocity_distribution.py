import mdconf as md
import matplotlib.pyplot as plt

fn_nw = '/home/michele/BeskowDiag/VelocityDistribution/NearWall'
fn_bk = '/home/michele/BeskowDiag/VelocityDistribution/ConfSnapshots'

v_bins, p_nw_z = md.velocity_distribution_binning( fn_nw, 200, 250, 7.5, 0.1, 'z' )
v_bins, p_nw_x = md.velocity_distribution_binning( fn_nw, 200, 250, 7.5, 0.1, 'x' )

v_bins, p_bk_z = md.velocity_distribution_binning( fn_bk, 200, 225, 7.5, 0.1, 'z' )
v_bins, p_bk_x = md.velocity_distribution_binning( fn_bk, 200, 225, 7.5, 0.1, 'x' )

plt.plot(v_bins, p_nw_z, 'r-', linewidth=2.0, label='z, wall')
plt.plot(v_bins, p_nw_x, 'b-', linewidth=2.0, label='x, wall')
plt.title('Dir. z vs x', fontsize=20.0)
plt.xlabel('v [nm/ps]', fontsize=20.0)
plt.ylabel('p [-1]', fontsize=20.0)
plt.legend(fontsize=20.0)
plt.xticks(fontsize=20.0)
plt.show()

plt.plot(v_bins, p_nw_z, 'r-', linewidth=2.0, label='z, wall')
plt.plot(v_bins, p_bk_z, 'b-', linewidth=2.0, label='z, bulk')
plt.title('Wall vs bulk', fontsize=20.0)
plt.xlabel('v [nm/ps]', fontsize=20.0)
plt.ylabel('p [-1]', fontsize=20.0)
plt.legend(fontsize=20.0)
plt.xticks(fontsize=20.0)
plt.show()
