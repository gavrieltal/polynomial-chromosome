class Expr # immutable Chromosome

  # an Expr is a polynomial represented by coefficients starting from x^0
  # e.g. [1, -2, 3] => 1 - 2x + 3x^2
  
  attr_reader :terms
  
  def initialize args
    if args.is_a?(Array)
      @terms = args
    elsif args.is_a?(Hash) && args.has_key?(:deg) && args.has_key?(:coeff_max)
      @terms = Array.new(args[:deg] + 1) {
        rand(args[:coeff_max]) * random_sign }
    else
      raise "args is an array of coeffs or a hash with :deg and :coeff_max"
    end
  end
      
  def mate_with expr, options = {}
    mutate_prob = options[:mutate_p] || 0.02
    abs_max     = options[:abs_max]  ||   -1
    degree      = [@terms.length, expr.terms.length].max
    crossover   = rand(degree)
    a = @terms[0...crossover] + expr.terms[crossover...degree]
    b = expr.terms[0...crossover] + @terms[crossover...degree]

    if abs_max.negative?
      abs_max = (a + b).map {|c| c.abs}.max - 1
    end
    
    if mutate_prob.positive?
      f = lambda {|i|
        if rand < mutate_prob then random_sign*rand(abs_max) else i end}
      a = a.map  {|i| f.call(i)}
      b = b.map  {|i| f.call(i)}
    end
    
    return [Expr.new(a), Expr.new(b)]
  end

  def to_s
    terms.each_with_index
      .map  {|v, d| v.to_s + "*x**" + d.to_s}
      .join ("+")
  end

  private

  def random_sign() rand(2).zero? ? 1 : -1 end
  
end


class Population
  
  attr_accessor :exprs, :inputs, :outputs, :cost_func, :bound_x, :bound_y,
                :mutate_p
  attr_reader   :scores, :graphed, :frame_num, :id, :abs_max
  
  def initialize(arr_inputs, arr_outputs, options = {})
    # available options: :num_exprs, :max_degree, :cost_function, :mutate_p
    num_exprs  = options[:num_exprs ] || 10
    max_degree = options[:max_degree] ||  3
    arr_inputs_max  = arr_inputs.max
    arr_outputs_max = arr_outputs.max

    if (arr_inputs.length != arr_outputs.length)
      raise "arr_inputs and arr_outputs must be of same length"
    end

    @abs_max     = [arr_inputs_max, arr_outputs_max].max + 1
    @exprs       = Array.new(num_exprs) {
      Expr.new({deg: max_degree, coeff_max: abs_max}) }
    @inputs      = arr_inputs
    @outputs     = arr_outputs
    @graphed     = false
    @frame_num   = 0
    @cost_func   = options[:cost_function] || (lambda {|x, y| (x - y)**2})
    @bound_x     = arr_inputs_max  + 1
    @bound_y     = arr_outputs_max + 1
    @id          = self.hash.abs
    @mutate_p    = options[:mutate_p] || 0.10  #probability of mutation
  end

  def inputs=  arr_inputs
    if arr_inputs.length != @outputs.length
      raise "arr_inputs must be same length as current outputs. Try set_values"
    else @inputs = arr_inputs
    end
  end

  def outputs= arr_outputs
    if arr_outputs.length != @inputs.length
      raise "arr_outputs must be same length as current inputs. Try set_values"
    else @outputs = arr_outputs
    end
  end

  def set_values arr_inputs, arr_outputs
    if arr_inputs.length != arr_outputs.length
      raise "arr_outputs must be same length as arr_inputs"
    else
      @inputs  = @arr_inputs
      @outputs = @arr_outputs
    end
  end

  def push_values input, output
    @inputs .push(input)
    @outputs.push(output)
    :ok
  end

  def push_expr arg = {deg: @exprs[0].terms.length, coeff_max: @abs_max}
    if arg.instance_of? Expr
      @exprs.push(arg)
    elsif arg.instance_of? Hash
      @exprs.push(Expr.new(arg))
    end
  end
  
  def scores
    if @scores.nil?
      @scores = @exprs.map {
        |e| (0...@inputs.length).map {
          |i| @cost_func.call(
            @outputs[i], e.terms.map.with_index { |coeff, exp|
              coeff*(@inputs[i]**exp) }.reduce(:+) ) }.reduce(:+) }
    else
      @scores
    end
  end
  
  def graph_exprs(savefile = "exprs.txt")
    if @graphed
      puts("this generation has already been graphed")
    else 
      File.open savefile, 'w' do |file|
        @exprs.map {|e| file.write(e.to_s + "\n")}
        file.close
      end

      command = "python graph_script.py #{savefile} #{id}" + 
                " #{frame_num} #{bound_x} #{bound_y}"
      puts(command)
      %x{#{command}}
      
      @graphed = true
    end
  end
  
  def mate_exprs # crossover_method
    # beforehand, calculate scores
    scores()
    
    # first, sort exprs by score
    @exprs = score_pairs().map { |pair| @exprs[pair[:index]] }
 
    # second, take exprs[2...-2] and mate them
    save_first = save_last = 2
    
    middle_exprs   = @exprs[save_first...-save_last]
    middle_indexes = (0...middle_exprs.length).to_a

    while !middle_exprs.empty?
      sample = middle_exprs.sample(2)
      index  = @exprs.length - middle_exprs.length - save_last

      if sample.length >= 2
        @exprs[index], @exprs[index + 1] =
                       sample[0].mate_with(sample[1],
                                           {mutate_prob: @mutate_p,
                                            abs_max: -1})
      else
        @exprs[index] = sample[0]
      end

      middle_exprs = middle_exprs - sample
    end

    # third, take exprs[0, 2], mate them, and replace worst scorers with them
    for i in (0...save_first).step(2)
      @exprs[-i-1], @exprs[-i-2] =
                    @exprs[i].mate_with(@exprs[i + 1],
                                        {mutate_prob: @mutate_p,
                                         abs_max: -1})
    end

    # fourth, implicitly do nothing to save_first's exprs

    # lastly, set Population's metadata to next generation/frame
    next_generation()
  end

  def best_fit # expr with lowest score
    best_pair = score_pairs()[0]
    @exprs[best_pair[:index]].to_s
  end
    
  def animate options = {}
    frames_per_sec = options[:frames_per_sec] ||  4
    fidelity       = options[:fildelity]      || 12
    rm             = if options[:rm].nil? then true else options[:rm] end
    
    %x{ffmpeg -r #{frames_per_sec} -i #{@id}_%d.png -pix_fmt yuv420p -r \
    #{fidelity} #{@id}.mp4}
    
    if rm
      %x{rm #{@id}*.png}
    end
    
    puts("your video is saved as ./#{@id}.mp4!")
  end

  def quick_movie frames = (@inputs.length**2).ceil
    puts("frames to produce: " + frames.to_s)
    frames.times {self.graph_exprs; self.mate_exprs}
    self.animate
    puts("calculated best fit expression: ")
    self.best_fit
  end
  
  def to_s
    print "inputs:  ", @inputs , "\n"
    print "outputs: ", @outputs, "\n"
    print "exprs:", "\n"
    @exprs.map { |e| e.to_s }
  end
    
  private

  def score_pairs
    scores()
    @scores
      .map.with_index { |score, i| { score: score, index: i } }
      .sort_by        { |pair| pair[:score]} 
  end
    
  def next_generation
    @graphed    = false
    @frame_num += 1
    @scores     = nil
  end
  
end
