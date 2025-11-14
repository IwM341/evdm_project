#pragma once

#include <variant>
#include <vector>
#include <memory>
#include <cmath>
#include <array>
#include <ranges>

namespace evdm {
	namespace detail{
		typedef uint_least64_t u64_t;
		u64_t morton_encode(u64_t x, u64_t y) {
			// Магические числа для 64-битного интерливинга
			x = (x | (x << 16)) & 0x0000FFFF0000FFFF;
			x = (x | (x << 8))  & 0x00FF00FF00FF00FF;
			x = (x | (x << 4))  & 0x0F0F0F0F0F0F0F0F;
			x = (x | (x << 2))  & 0x3333333333333333;
			x = (x | (x << 1))  & 0x5555555555555555;
			
			y = (y | (y << 16)) & 0x0000FFFF0000FFFF;
			y = (y | (y << 8))  & 0x00FF00FF00FF00FF;
			y = (y | (y << 4))  & 0x0F0F0F0F0F0F0F0F;
			y = (y | (y << 2))  & 0x3333333333333333;
			y = (y | (y << 1))  & 0x5555555555555555;
			
			return x << 1 | y;
		}
		std::pair<u64_t, u64_t> morton_decode(u64_t z) {
			u64_t x = z & 0x5555555555555555;
			u64_t y = (z >> 1) & 0x5555555555555555;
			
			x = (x | (x >> 1))  & 0x3333333333333333;
			x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0F;
			x = (x | (x >> 4))  & 0x00FF00FF00FF00FF;
			x = (x | (x >> 8))  & 0x0000FFFF0000FFFF;
			x = (x | (x >> 16)) & 0x00000000FFFFFFFF;
			
			y = (y | (y >> 1))  & 0x3333333333333333;
			y = (y | (y >> 2))  & 0x0F0F0F0F0F0F0F0F;
			y = (y | (y >> 4))  & 0x00FF00FF00FF00FF;
			y = (y | (y >> 8))  & 0x0000FFFF0000FFFF;
			y = (y | (y >> 16)) & 0x00000000FFFFFFFF;
			
			return {y, x};
		}
	};

	struct Nothing_t{
		template <typename Serializer_t>
		auto Serialize(Serializer_t && ser) const{
			return grob::Serialize(false,ser);
		}
		template <typename Object_t,typename Deserializer_t>
		static Nothing_t DeSerialize(Object_t const& obj,Deserializer_t && deser){
			return Nothing_t{};
		}
	};

	template <typename Fp_t,typename T,typename U>
	class QuadtreeNode{
	private:
		typedef uint_least64_t u64_t;
		//(0,0)=(left down) (0,1)=(left up), (1,0) = (right down), (1,1) = (right up)
		std::array<std::unique_ptr<QuadtreeNode>, 4> children = { nullptr,nullptr,nullptr,nullptr }; // NW, NE, SW, SE
		QuadtreeNode* _parent = nullptr;
		u64_t _level = 0;
		u64_t cx0, cx1, cy0, cy1;
		Fp_t _x0, _x1, _y0, _y1;
		U m_value;
		std::array<std::shared_ptr<T>, 4> corner_values; // NW, NE, SW, SE corners
		bool _is_leaf;

		template <typename Fp1_t,typename T1,typename U1,size_t maxlevel>
		friend struct Quadtree;
		QuadtreeNode(
			u64_t cx0, u64_t cx1, u64_t cy0, u64_t cy1,
			Fp_t _x0, Fp_t _x1, Fp_t , Fp_t _y1,
			QuadtreeNode* m_parent = nullptr, u64_t m_level = 0)
			: cx0(cx0), cx1(cx1), cy0(cy0), cy1(cy1),
			_x0(_x0), _x1(_x1),
			_y0(_y0), _y1(_y1),
			_parent(m_parent), _level(m_level), _is_leaf(true) {}


	public:

		inline U & uvalue(){
			return m_value;
		}
		inline U const& uvalue()const {
			return m_value;
		}

#define DECLARE_CONST_PROPERTY(Type,prop) \
	inline Type prop() const {return _##prop;}
		DECLARE_CONST_PROPERTY(Fp_t, x0)
			DECLARE_CONST_PROPERTY(Fp_t, x1)
			DECLARE_CONST_PROPERTY(Fp_t, y0)
			DECLARE_CONST_PROPERTY(Fp_t, y1)
#undef DECLARE_CONST_PROPERTY
		inline grob::Point<grob::Rect<Fp_t>,grob::Rect<Fp_t>> 
			rect() const{
				return { {_x0,_x1},{_y0,_y1}};
			}
		inline size_t level()const {
			return _level;
		}
		inline bool is_leaf() const{
			return _is_leaf;
		}
		inline QuadtreeNode& child(size_t i) {
			if (_is_leaf) {
				throw std::runtime_error("leaf has no child");
			}
			return *children[i];
		}
		inline QuadtreeNode const& child(size_t i)const {
			if (_is_leaf) {
				throw std::runtime_error("leaf has no child");
			}
			return *children[i];
		}

		inline QuadtreeNode& parent() {
			if (_parent == nullptr) {
				throw std::runtime_error("root has no parent");
			}
			return *_parent;
		}
		inline QuadtreeNode* parentp() {
			return _parent;
		}
		inline const QuadtreeNode* parentp()const {
			return _parent;
		}
		inline bool contain(Fp_t x, Fp_t y)const {
			return _x0 <= x && x <= _x1 && _y0 <= y && y <= _y1;
		}
		inline std::pair<Fp_t, Fp_t> coeffs(Fp_t x, Fp_t y)const {
			return { (x1() - x) / (x1() - x0()),(y1() - y) / (y1() - y0()) };
		}
		inline std::tuple<const T&,const T&,const T&,const T&> values()const {
			return {*corner_values[0],*corner_values[1],*corner_values[2],*corner_values[3]};
		}
		template <typename U1>
		using tuple4 = std::tuple<U1,U1,U1,U1>; 

		template <typename ForEachFunc_t>
		inline tuple4<std::invoke_result_t<ForEachFunc_t,const T&>> 
			values_view(ForEachFunc_t && F)const {
			return {
				F(*corner_values[0]),
				F(*corner_values[1]),
				F(*corner_values[2]),
				F(*corner_values[3])
			};
		}
		inline const T& operator ()(size_t i) const {
			return *corner_values[i];
		}
		inline T& operator ()(size_t i) {
			return *corner_values[i];
		}
		inline const T& operator ()(size_t i, size_t j) const {
			return *corner_values[2 * i + j];
		}
		inline T& operator ()(size_t i, size_t j) {
			return *corner_values[2 * i + j];
		}
		std::array<u64_t,4> coords()const{
			return {cx0,cx1,cy0,cy1};
		}
	};


	template<typename Fp_t, typename T,typename U = Nothing_t, size_t max_lvl = 30>
	struct Quadtree {
		typedef uint_least64_t u64_t;
		using value_type = T;
		using coord_type =  Fp_t;
		using Node = QuadtreeNode<Fp_t,T,U>;


		std::unique_ptr<Node> m_root;
		typedef std::map<u64_t, std::shared_ptr<T>> PointValues_t;
		PointValues_t  point_registry;

		Fp_t X0, X1, Y0, Y1;
		//std::atomic<size_t> max_level_dyn = 0;
		size_t leaf_count = 1;
		size_t node_count = 1;

		constexpr static u64_t _1L = ((u64_t)1);
		constexpr static u64_t _3L = ((u64_t)3);
		constexpr static u64_t max_c = _1L << max_lvl;
		constexpr static u64_t max_c_inner = (_1L << max_lvl) - 1;

		Fp_t dX_inv() const {
			return 1 / (X1 - X0);
		}
		Fp_t dY_inv() const {
			return 1 / (Y1 - Y0);
		}
		// Хеш функция для точки
		u64_t point_coord(Fp_t x, Fp_t y) const {
			u64_t hx = x_coord(x);
			u64_t hy = y_coord(y);
			return (hx << (max_lvl + 1)) | hy;
		}

		u64_t x_coord(Fp_t x)const {
			auto max_c_d = max_c;
			return static_cast<u64_t>(((x - X0) * dX_inv()) * max_c);
		}
		u64_t y_coord(Fp_t y)const {
			auto max_c_d = max_c;
			return static_cast<uint_least64_t>(((y - Y0) * dY_inv()) * max_c);
		}

		Fp_t x_val(uint64_t xc)const {
			Fp_t a = (Fp_t(xc)) / max_c;
			return X0 + a * (X1 - X0);
		}
		Fp_t y_val(u64_t yc)const {
			Fp_t b = (Fp_t(yc)) / max_c;
			return Y0 + b * (Y1 - Y0);
		}
		inline static u64_t get_key(u64_t xc, u64_t yc) {
			return (xc << (max_lvl + 1)) | yc;
		}
		inline static std::pair<u64_t, u64_t> xcyc(u64_t key) {
			return { key >> (max_lvl + 1),key & max_c_inner };
		}

		// Получение или создание значения для точки
		template <typename Value_Gen_t>
		std::shared_ptr<T> get_or_create_value(
			u64_t xc, u64_t yc, 
			Value_Gen_t&& v_creator) {
			auto m_point = get_key(xc, yc);
			{
				
				if (auto it = point_registry.find(m_point); it != point_registry.end()) {
					return it->second;
				}
			}
			auto value = std::make_shared<T>(v_creator());
			{
				point_registry[m_point] = value;
			}

			return value;
		}

		

		// Наследование значений углов от родителя
		template <typename Fun_t>
		void inherit_corner_values(Node& node, Fun_t&& F) {

			node.corner_values[0] = get_or_create_value(node.cx0, node.cy0,
				[x = node.x0(), y = node.y0(), &F]() {return F(x, y); });

			node.corner_values[1] = get_or_create_value(node.cx0, node.cy1,
				[x = node.x0(), y = node.y1(), &F]() {return F(x, y); });

			node.corner_values[2] = get_or_create_value(node.cx1, node.cy0,
				[x = node.x1(), y = node.y0(), &F]() {return F(x, y); });

			node.corner_values[3] = get_or_create_value(node.cx1, node.cy1,
				[x = node.x1(), y = node.y1(), &F]() {return F(x, y); });

		}

		void inherit_existing_values(Node& node) {
			auto ExistingGetter = [cx = node.cx0,cy = node.cy0]()->T {
				std::ostringstream Err;
				Err << "corner value not found for node (";
				Err << cx << ", " << cy << ")";
				throw std::runtime_error("corner value not found");
			};
			node.corner_values[0] = get_or_create_value(
				node.cx0, node.cy0,ExistingGetter);

			node.corner_values[1] = get_or_create_value(
				node.cx0, node.cy1,ExistingGetter);

			node.corner_values[2] = get_or_create_value(
				node.cx1, node.cy0,ExistingGetter);

			node.corner_values[3] = get_or_create_value(
				node.cx1, node.cy1,ExistingGetter);

		}
		void inherit_existing_values_recursive(Node& node){
			inherit_existing_values(node);
			if(!node.is_leaf()){
				for(size_t i=0;i<4;++i){
					inherit_existing_values_recursive(node.child(i));
				}
			}
		}

		std::unique_ptr<Node> create_node(
			u64_t x0, u64_t x1, u64_t y0, u64_t y1,
			Node* parent = nullptr, int level = 0)
		{
			
			return std::make_unique<Node>(
				Node(
					x0, x1, y0, y1,
					x_val(x0),x_val(x1),y_val(y0),y_val(y1),
					parent, level));
		}

		const Quadtree& cthis()const {
			return *this;
		}

		private:
		Quadtree(Fp_t X0 = 0, Fp_t X1 = 1, Fp_t Y0 = 0, Fp_t Y1 = 1)
			:X0(X0), X1(X1), Y0(Y0), Y1(Y1)
		{
			static_assert(max_lvl <= 31, "max_lvl should be less 31");
			if (X1 <= X0 || Y1 <= Y0) {
				throw std::runtime_error("X1 <= X0 || Y1 <= Y0");
			}
			m_root = create_node(0, max_c, 0, max_c);
		}
		class StubFunc_t{};
	public:
        size_t leafs() const {
            return leaf_count;
        }
        size_t nodes() const {
            return node_count;
        }
		template <typename FillFunc,
			typename = std::enable_if_t<std::is_invocable_v<FillFunc,Fp_t,Fp_t>>
		>Quadtree(FillFunc&& F, Fp_t X0 = 0, Fp_t X1 = 1, Fp_t Y0 = 0, Fp_t Y1 = 1)
			:X0(X0), X1(X1), Y0(Y0), Y1(Y1)
		{
			static_assert(max_lvl <= 31, "max_lvl should be less 31");
			if (X1 <= X0 || Y1 <= Y0) {
				throw std::runtime_error("X1 <= X0 || Y1 <= Y0");
			}
			m_root = create_node(0, max_c, 0, max_c);
			inherit_corner_values(*m_root, F);
		}
		
		Node& root() {
			return *m_root;
		}
		Node const& root()const {
			return *m_root;
		}
		template <typename Value_Gen_t>
		void refine(Node& node, size_t num_levels, Value_Gen_t&& F){
			if (num_levels == 0) {
				return;
			}
			else {
				if (node.is_leaf()) {
					if constexpr(
						std::is_same_v<StubFunc_t,std::decay_t<Value_Gen_t>>
					){
						split_node_unsafe(node);
					} else{
						split_node(node, F);
					}
					for (auto& ch : node.children) {
						refine(*ch, num_levels - 1, F);
					}
				}
				else {
					for (auto& ch : node.children) {
						refine(*ch, num_levels, F);
					}
				}
			}
		}

		template <typename Value_Gen_t>
		void refine_to_level(Node& node, size_t max_level, Value_Gen_t&& F) {
			if(max_level > node.level()){
				refine(node,max_level-node.level(),F);
			}
		}

		void refine_unsave(Node& node, size_t num_levels){
			return refine(node,num_levels,StubFunc_t{});
		}
		void refine_to_level_unsave(Node& node, size_t num_levels){
			return refine_to_level(node,num_levels,StubFunc_t{});
		}
		/// @brief splits current node without putting values
		/// @param node 
		void split_node_unsafe(Node & node){
						if (!node.is_leaf()) return;
			if (node.level() >= max_lvl) {
				throw std::runtime_error("max level reached");
			}



			u64_t cx_n = (node.cx0 + node.cx1) / 2;
			u64_t cy_n = (node.cy0 + node.cy1) / 2;

			// Создаем дочерние узлы
			u64_t new_level = node.level() + 1;
			auto n00 = create_node(node.cx0, cx_n, node.cy0, cy_n, &node, new_level);
			auto n01 = create_node(node.cx0, cx_n, cy_n, node.cy1, &node, new_level);
			auto n10 = create_node(cx_n, node.cx1, node.cy0, cy_n, &node, new_level);
			auto n11 = create_node(cx_n, node.cx1, cy_n, node.cy1, &node, new_level);
			node.children[0] = std::move(n00);
			node.children[1] = std::move(n01);
			node.children[2] = std::move(n10);
			node.children[3] = std::move(n11);
			node._is_leaf = false;
			//max_level_dyn.compare_exchange_strong = std::max(max_level_dyn, new_level);
			leaf_count += (3);
			node_count += (4);
		}
		// Разделение узла
		template <typename Value_Gen_t>
		void split_node(Node& node, Value_Gen_t&& value_func) {
			split_node_unsafe(node);

			// Наследуем значения углов от родителя
			inherit_corner_values(*node.children[0], value_func);
			inherit_corner_values(*node.children[1], value_func);
			inherit_corner_values(*node.children[2], value_func);
			inherit_corner_values(*node.children[3], value_func);
		}

		//missing keys vector    
		//missing x,y vector
		std::vector<std::tuple<u64_t,Fp_t,Fp_t>>
			get_missing_values_unsafe() const
		{
			typedef std::tuple<u64_t,Fp_t,Fp_t> _Tp;
			struct m_compare{
				inline constexpr bool operator()(const _Tp& lhs, const _Tp& rhs)const{
					return std::get<0>(lhs) < std::get<0>(rhs); 
				}
			};
			std::set<std::tuple<u64_t,Fp_t,Fp_t>,m_compare> missing_keys;
			for(const Node & N : as_node_container()){
				for(size_t i=0;i<4;++i){
					auto cx = i/2 ? N.cx0 : N.cx1;
					auto cy = i%2 ? N.cy0 : N.cy1;
					
					auto m_key = get_key(cx,cy);
					if(
						point_registry.find(m_key) == 
						point_registry.end()
					){
						auto xp = i/2 ? N.x0() : N.x1();
						auto yp = i%2 ? N.y0() : N.y1();
					
						missing_keys.emplace(m_key,xp,yp);
					}
				}
			}
			return {missing_keys.begin(),missing_keys.end()};
		}
		void insert_missing_values(
			std::vector<
				std::pair<u64_t,std::shared_ptr<T>>
			> m_values)
		{
			point_registry.merge(PointValues_t(m_values.begin(),m_values.end()));
			inherit_existing_values_recursive(root());
		}

		/// @brief refine nodes
		/// @param r_cond if r_cond(node) == true node will be refined
		/// @param value_provider gives from vector<pair(x_i,y_i)> vector<std::shared_ptr<T>>
		/// @return number of new nodes
		template <typename RefineCondition_t, typename XYFuncVectorization_t>
		size_t refine_vectorized(RefineCondition_t r_cond, XYFuncVectorization_t && value_provider,
			size_t m_level = 1,size_t max_new_nodes = std::numeric_limits<size_t>::max()){
			size_t new_nodes = 0;
			for (Node & N : *this){
				size_t delta = (1<<m_level) - 1;
				if(max_new_nodes >= delta && r_cond(N)){
					refine_unsave(N,m_level);
					max_new_nodes -=delta;
					new_nodes += delta;
				}
			}
			auto missings_tp = get_missing_values_unsafe();
			auto as_pair = missings_tp | std::ranges::views::transform(
				[](auto const & tp) -> std::pair<Fp_t,Fp_t> {return {std::get<1>(tp),std::get<2>(tp)};}
			);
			std::vector<std::pair<Fp_t,Fp_t>> missings(as_pair.begin(),as_pair.end());
			std::vector<std::shared_ptr<T>> Values = value_provider(missings);
			if( missings.size() != missings.size()) {
				throw std::runtime_error("Quadtree::refine_vectorized : missings.size() != missings.size()");
			}
			std::vector< std::pair<u64_t,std::shared_ptr<T>> > m_values(missings.size());
			for(size_t i=0;i<m_values.size();++i){
				m_values[i] = {std::get<0>(missings_tp[i]),std::move(Values[i])};
			}
			insert_missing_values(m_values);
			return new_nodes;
		}

		// Быстрый поиск узла содержащего точку
		inline Node const& find_node_coord(u64_t cx, u64_t cy,Node & _tmp_node) const {
			Node * tmp_node = &_tmp_node;
			while (!tmp_node->is_leaf()) {
				u64_t ix = 1 & (cx >> (max_lvl - tmp_node->_level - 1));
				u64_t iy = 1 & (cy >> (max_lvl - tmp_node->_level - 1));
				tmp_node = &tmp_node->child(ix << 1 | iy);
			}
			return *tmp_node;
		}
		// Быстрый поиск узла содержащего точку
		inline Node & find_node_coord(u64_t cx, u64_t cy,Node & tmp_node)  {
			return const_cast<Node &>( (cthis().find_node_coord(cx,cy,tmp_node)));
		}

		// Быстрый поиск узла содержащего точку
		inline Node const& find_node_coord(u64_t cx, u64_t cy) const {
			return find_node_coord(cx,cy,*m_root.get());
		}
		inline Node& find_node_coord(u64_t cx, u64_t cy) {
			return const_cast<Node&>(cthis().find_node_coord(cx, cy));
		}

		Node const& find_node(Fp_t x, Fp_t y) const {
			x = std::clamp(x, X0, X1);
			y = std::clamp(y, Y0, Y1);

			u64_t cx = std::min(x_coord(x), max_c_inner);
			u64_t cy = std::min(y_coord(y), max_c_inner);
			return find_node_coord(cx, cy);
		}
		inline Node& find_node(Fp_t x, Fp_t y) {
			return const_cast<Node&>(cthis().find_node(x, y));
		}

		bool contain(Fp_t x, Fp_t y) const {
			return root().contain(x, y);
		}

		template <typename Node_t>
		struct node_iterator_t {
			std::stack<Node_t*> m_stack;
			
			node_iterator_t operator++(){
				if(m_stack.empty()) return *this;
				Node_t* m_current = m_stack.top();
				m_stack.pop();
				if( !m_current->is_leaf()){
					for (size_t i=0;i<4;++i){
						m_stack.push(& (m_current->child(i)));
					}
				}
				return *this;
			}
			Node_t & operator*() { return *m_stack.top(); }
			Node_t* operator->() { return m_stack.top(); }
			operator bool() const{
				return !m_stack.empty();
			}
			inline bool operator==(const node_iterator_t& other) const {
				return m_stack.empty();
			}
			inline bool operator!=(const node_iterator_t& other) const {
				return !(*this == other);
			}
			node_iterator_t(Node_t* _root) {
				if (_root) m_stack.push(_root);
			}
			node_iterator_t(){}
		};

		template <typename Quadtree_t>
		struct node_container{
			Quadtree_t & m_this;
			typedef std::remove_reference_t<
				decltype(m_this.root())
			> Node_t;

			node_iterator_t<Node_t> begin(){
				return {&m_this.root()};
			}
			node_iterator_t<Node_t> end(){
				return {&m_this.root()};
			}
			node_iterator_t<const Node_t> begin()const{
				return {&m_this.root()};
			}
			node_iterator_t<const Node_t> end()const{
				return {&m_this.root()};
			}

			node_iterator_t<const Node_t> cbegin()const{
				return {&m_this.root()};
			}
			node_iterator_t<const Node_t> cend()const{
				return {&m_this.root()};
			}
			
		};

		node_container<Quadtree> as_node_container(){
			return {*this};
		}
		node_container<const Quadtree> as_node_container()const{
			return {*this};
		}
		
		template <typename Node_t>
		struct leaf_iterator_t {
			Node_t* m_node;

			leaf_iterator_t(Node_t* m_node) :m_node(m_node) {}
			operator bool() const{
				return m_node != nullptr;
			}
			inline bool operator ==(leaf_iterator_t const& other)const {
				return m_node == other.m_node;
			}
			inline leaf_iterator_t& operator++() {
				u64_t morton_code = detail::morton_encode(
					m_node->cx0, m_node->cy0
				);
				/*
				println("tmpcode");
				println(std::bitset<8>(morton_code));
				println(std::bitset<4>(m_node->cx0));
				println(std::bitset<4>(m_node->cy0));
				println("tmp codes");
				std::cout << std::bitset<31>(m_node->cx0) <<
					"\n" << std::bitset<31>(m_node->cy0) << std::endl;
				std::cout << std::bitset<62>(morton_code)<< std::endl;

				u64_t morton_code_new = detail::morton_encode(
							m_node->_parent->children[1]->cx0,
							m_node->_parent->children[1]->cy0
						);
				println("desire code");
				std::cout << std::bitset<31>(
					m_node->_parent->children[1]->cx0) << std::endl;

				std::cout << std::bitset<31>(
					m_node->_parent->children[1]->cy0) << std::endl;

				std::cout << std::bitset<62>(morton_code_new)<< std::endl;
				*/
				u64_t new_code = morton_code + (_1L << (2 * (max_lvl - m_node->_level)));

				//println("new code");
				//println(std::bitset<8>(new_code));
				/*println("new code");
				std::cout << std::bitset<62>(new_code)<< std::endl;
				*/
				//assert(cy == );
				Node_t* common_parent = find_common_parent(
					m_node, morton_code, new_code,
					(_3L) << (2 * (max_lvl - m_node->_level))
				);
				/*println("mask");
				std::cout << std::bitset<62>(
					(_3L) << (2*(max_lvl - m_node->_level))
				)<< std::endl;*/
				if (common_parent == nullptr) {
					m_node = nullptr;
					return *this;
				}

				m_node = find_leaf(common_parent, new_code,
					2 * (max_lvl - common_parent->_level - 1));

				//assert(m_node == &find_node_coord(cx,cy));
				return *this;
			};
			inline Node_t& operator *() {
				return *m_node;
			}
			inline Node_t* operator ->() {
				return m_node;
			}

		private:
			friend struct Quadtree;
			static inline Node_t* find_common_parent(
				Node_t* tmp_node, u64_t zcode, u64_t tmp_code, u64_t level_mask) {
				while ((zcode & level_mask) != (tmp_code & level_mask)) {
					tmp_node = tmp_node->_parent;
					level_mask = level_mask << 2;
				}
				return tmp_node;
			}
			static inline Node_t* find_leaf(
				Node_t* parent, u64_t m_code, u64_t m_level) {
				while (!parent->is_leaf()) {
					size_t index = (m_code >> m_level) & _3L;
					parent = &(parent->child(index));
					m_level -= 2;
				}
				return parent;
			}
		};

		leaf_iterator_t<Node> begin() {
			return { &find_node_coord(0,0) };
		}

		leaf_iterator_t<Node> end() {
			return { nullptr };
		}
		leaf_iterator_t<const Node> cbegin() const{
			return { &find_node_coord(0,0) };
		}

		leaf_iterator_t<const Node> cend() const{
			return { nullptr };
		}
		leaf_iterator_t<const Node> begin() const{
			return { &find_node_coord(0,0) };
		}

		leaf_iterator_t<const Node> end() const{
			return { nullptr };
		}

		template <typename Serializator_t>
		auto Serialize(Serializator_t && S)const{
			typedef std::tuple<size_t,u64_t,u64_t> NodeCoord;

			//                 level  xcoord ycoord
			std::vector<NodeCoord> LeafCoords;
			LeafCoords.reserve(leaf_count);
			std::vector<std::reference_wrapper<const U>> ValuesU;
			ValuesU.reserve(node_count);

			for( const Node & N: as_node_container() ){
				ValuesU.push_back( N.uvalue());
			}
			
			for (const Node & N: *this ){
				LeafCoords.push_back({N._level,N.cx0,N.cy0});
			}


			std::vector<
				std::pair<
					u64_t, // key
					std::reference_wrapper<const T> // value
				>
			> Values;

			Values.reserve(point_registry.size());


			for(auto &[key,value] : point_registry){
				Values.push_back({key,*value});
			}
			
			//static_assert(false);
			return S.MakeDict(
				grob::make_vector("Box","Leafs","Values","ValuesU"),
				grob::make_vector(
					grob::Serialize(std::make_tuple(X0,X1,Y0,Y1),S),
					grob::Serialize(std::move(LeafCoords),S),
					grob::Serialize(std::move(Values),S),
					grob::Serialize(std::move(ValuesU),S)
				)
			);
		}

		template <typename Object_t,typename Dser_t>
		static Quadtree DeSerialize(Object_t && Obj,Dser_t && DS){
			typedef std::tuple<size_t,u64_t,u64_t> NodeCoord;
			typedef std::vector<NodeCoord> NCV_t;
			typedef std::pair<u64_t,T> KeyVal_t;
			typedef std::vector<KeyVal_t> MapVecVal_t;
			typedef std::tuple<Fp_t,Fp_t,Fp_t,Fp_t,NCV_t,MapVecVal_t> DS_Result_t;

			auto Properties = grob::GetProperties(
				grob::make_vector("Box","Leafs","Values","ValuesU"),DS,Obj
			);
			auto m_box = grob::DeSerialize<std::tuple<Fp_t,Fp_t,Fp_t,Fp_t>>(
				Properties[0],DS
			);
			Fp_t m_X0 =  std::get<0>(m_box);
			Fp_t m_X1 =  std::get<1>(m_box);
			Fp_t m_Y0 =  std::get<2>(m_box);
			Fp_t m_Y1 =  std::get<3>(m_box);
			
			NCV_t Nodes = grob::DeSerialize<NCV_t>(Properties[1],DS);
			MapVecVal_t m_map_vec_pre = grob::DeSerialize<MapVecVal_t>(Properties[2],DS);

			std::vector< std::pair<u64_t,std::shared_ptr<T>>> m_map_vec;
			m_map_vec.reserve(m_map_vec_pre.size());
			for(auto && v : m_map_vec_pre){
				m_map_vec.push_back({ v.first,std::make_shared<T>(std::move(v.second)) });
			}

			std::vector<U> m_values_u = 
				grob::DeSerialize<std::vector<U>>(Properties[3],DS);

			Quadtree m_tree(m_X0,m_X1,m_Y0,m_Y1);
			m_tree.point_registry = 
				std::map<u64_t, std::shared_ptr<T>>( m_map_vec.begin(),m_map_vec.end());
			
			for (auto [level,cx,cy] : Nodes){
				Node * best_node = &m_tree.find_node_coord(cx,cy);
				while(best_node->level() < level){
					m_tree.split_node_unsafe(*best_node);
					best_node = &m_tree.find_node_coord(cx,cy,*best_node);
				}
			}
			//Quadtree::Quadtree;
			m_tree.inherit_existing_values_recursive(m_tree.root());

			size_t i = 0;
			for (Node & N : m_tree.as_node_container()){
				N.uvalue() = std::move( m_values_u[i] );
				++i;
			}

			return m_tree;
		}
	};
	

};