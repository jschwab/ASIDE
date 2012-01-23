# plot_aside.rb

require 'Tioga/FigureMaker'
require './plot_styles.rb'

class MyPlots
  include Math
  include Tioga
  include FigureConstants
  include MyPlotStyles
  
  def t
    @figure_maker
  end

  def initialize
    @figure_maker = FigureMaker.default
    
    t.save_dir = 'plots_out'
    @data_filename = "../exe/aside.out"
    @swifter_filename = "/Users/jschwab/Classes/Fall_2011/AY250/integrators/swifter/bin/follow.out"

    # configure data array
    @time = Dvector.new
    @a = Dvector.new
    @e = Dvector.new
    @i = Dvector.new
    @o = Dvector.new
    @w = Dvector.new
    @f = Dvector.new

    @id = Dvector.new
    @apo = Dvector.new
    @peri = Dvector.new

    @swifter_time = Dvector.new
    @swifter_a = Dvector.new
    @swifter_e = Dvector.new
    @swifter_i = Dvector.new
    @swifter_o = Dvector.new
    @swifter_w = Dvector.new
    @swifter_f = Dvector.new
    @swifter_id = Dvector.new
    @swifter_apo = Dvector.new
    @swifter_peri = Dvector.new


    # conversion factor
    @RAD2DEG = 57.29577951308232087679815

    @data_array = [@time,  @a,@e,@i,@o,@w,@f]
    @swifter_array = [@swifter_time, @swifter_id, 
                      @swifter_a,@swifter_e,@swifter_i,
                      @swifter_o,@swifter_w,@swifter_f,\
                      @swifter_apo,@swifter_peri]

    # props
    @margin = 0.05

    @compare = TRUE

    # define figures
    t.def_figure("semimajoraxis") { semimajoraxis } 
    t.def_figure("eccentricity") { eccentricity } 
    t.def_figure("inclination") { inclination } 
    t.def_figure("ascendingnode") { ascendingnode } 
    t.def_figure("argumentofpericenter") { argumentofpericenter } 

    t.def_figure("orbel") {orbel}
    t.def_enter_page_function { enter_page }

  end

  def enter_page
    set_default_plot_style
    t.default_enter_page_function
  end


  def plot_boundaries(xs,ys,margin,
                      xmin=nil,xmax=nil,ymin=nil,ymax=nil,
                      reverse_xaxis=false,reverse_yaxis=false)
    xmin = xs.min if xmin == nil
    xmax = xs.max if xmax == nil
    ymin = ys.min if ymin == nil
    ymax = ys.max if ymax == nil
    width = (xmax == xmin)? 1 : xmax - xmin
    height = (ymax == ymin)? 1 : ymax - ymin
    left_boundary = xmin - margin * width
    right_boundary = xmax + margin * width
    top_boundary = ymax + margin * height
    bottom_boundary = ymin - margin * height
    if reverse_xaxis
      tmp = left_boundary; left_boundary = right_boundary; right_boundary = tmp
    end
    if reverse_yaxis
      tmp = top_boundary; top_boundary = bottom_boundary; bottom_boundary = tmp
    end
    return [ left_boundary, right_boundary, top_boundary, bottom_boundary ]
  end


  def read_data

    if not @read_once
      read_aside
      if @compare
        read_swifter
      end
    @read_once = TRUE
    end

  end

  def read_aside
    puts @data_filename
    Dvector.read(@data_filename, @data_array)
    @tmin = @time.min
    @tmax = @time.max
  end

  def read_swifter
    Dvector.read(@swifter_filename, @swifter_array)
  end

  def make_panel(ys, ylabel, ycompare = nil)

    read_data

    xs = @time
    xlabel = "Time [s]"

    t.set_aspect_ratio(0.666666)    

    t.show_plot(plot_boundaries(xs,ys,@margin)) {

      t.show_xlabel(xlabel)
      t.show_ylabel(ylabel)      

      t.show_polyline(xs,ys)
      if @compare and ycompare != nil
        t.show_polyline(@swifter_time, ycompare, Red, nil, LINE_TYPE_DASH)
      end
    }

  end

  def semimajoraxis
    make_panel(@a,"a",@swifter_a)
  end


  def eccentricity
    make_panel(@e,"e",@swifter_e)
  end

  def inclination
    ys = @i* @RAD2DEG
    make_panel(ys,"$i$ [deg]",@swifter_i)
  end

  def ascendingnode
    ys = @o * @RAD2DEG
    make_panel(ys,'$\Omega$ [deg]', @swifter_o)
  end

  def argumentofpericenter
    ys = @w * @RAD2DEG
    make_panel(ys,'$\omega$ [deg]', @swifter_w)
  end

  def jacobiconstant
    ys = 0.5 / @a + (@a * (1-@e*@e)).sqrt()
    make_panel(ys, "Jacobi Constant")
  end


  def orbel

    t.rescale(0.67)
    column_margin = 0.1
    t.subplot(t.column_margins(
        'num_columns' => 3, 'column' => 1,
        'column_margin' => column_margin)) {
      t.subplot(t.row_margins(
          'num_rows' => 2, 'row' => 1,
          'row_margin' => column_margin)) { semimajoraxis }
      t.subplot(t.row_margins(
          'num_rows' => 2, 'row' => 2,
          'row_margin' => column_margin)) { ascendingnode }
    }

    t.subplot(t.column_margins(
        'num_columns' => 3, 'column' => 2,
        'column_margin' => column_margin)) {
      t.subplot(t.row_margins(
          'num_rows' => 2, 'row' => 1,
          'row_margin' => column_margin)) { eccentricity }
      t.subplot(t.row_margins(
          'num_rows' => 2, 'row' => 2,
          'row_margin' => column_margin)) { argumentofpericenter }
    }

    t.subplot(t.column_margins(
        'num_columns' => 3, 'column' => 3,
        'column_margin' => column_margin)) {
      t.subplot(t.row_margins(
          'num_rows' => 2, 'row' => 1,
          'row_margin' => column_margin)) { inclination }
      t.subplot(t.row_margins(
          'num_rows' => 2, 'row' => 2,
          'row_margin' => column_margin)) { jacobiconstant }
    }

  end

end

MyPlots.new
